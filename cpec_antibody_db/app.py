"""CPEC antibody database: guided interface.

Run with:  streamlit run app.py   (or double-click run.bat on Windows, ./run.sh on Mac or Linux)
"""
import io
import urllib.parse
import zipfile

import pandas as pd
import streamlit as st

import autofill
import db
import labels as L
import panels
import ranking

st.set_page_config(page_title="CPEC Antibody Database", layout="wide")
db.ensure_db()

CELLS = "stem-cell-derived CPECs"
TAG = {"ok": "[OK]", "check": "[Check]", "problem": "[Problem]"}
ROLE_OPTIONS = list(L.ROLE.values())
ROLE_CODE = {v: k for k, v in L.ROLE.items()}
HELP = {
    "Score": "Higher means a better marker to test first. The Score Breakdown on each marker shows the parts.",
    "Role": "Keep marks CPECs you want to collect. Remove marks cells you want to discard, such as leftover stem cells.",
    "Estimated Live-Cell Access": "Estimated from UniProt topology: how much of the protein sits outside the cell. "
                                  "Good means a GPI anchor or a large outside segment.",
    "Longest Outside Segment (aa)": "Length in amino acids of the longest stretch outside the cell. "
                                    "More length gives an antibody more to bind.",
    "Specificity": "Whether other cell types also make this protein. Shared or Found Everywhere means it needs a partner marker.",
    "Choroid Plexus RNA (nTPM)": "RNA level in choroid plexus from the Human Protein Atlas. RNA does not prove surface protein.",
    "Papers Found": "Papers on Europe PMC that mention the marker with choroid plexus and sorting or surface terms.",
    "Antibodies on File": "Antibodies saved in this database for the marker, not counting ones aimed inside the cell.",
    "Evidence on File": "Only what is recorded: your tests, live-cell sorts reported elsewhere, or the basis for the CPEC link.",
    "Immunogen Location": "Where the antibody's immunogen lies compared with the protein topology. It must be outside the cell.",
    "Problems": "Datasheet checks that make the antibody unlikely to work for live-cell sorting.",
    "To Check": "Datasheet items that are unknown or need confirming with the vendor.",
    "RRID": "Research Resource Identifier: a permanent ID for the antibody product that you can cite in papers.",
    "Catalog": "The vendor's catalog number.", "Citations": "Papers that cite the product. Use shows popularity, not live-cell success.",
    "Status": "Clone Name Only, Details Missing, Vendor Listed, or Passed in Your Lab.",
}


# ---------- helpers ----------
def flash(msg):
    st.session_state["flash"] = msg
    st.rerun()


def table(df, key=None):
    cfg = {c: st.column_config.Column(help=HELP[c]) for c in df.columns if c in HELP}
    st.dataframe(df, hide_index=True, width="stretch", column_config=cfg, key=key)


def num(series):
    return pd.to_numeric(series, errors="coerce")


def vendor_links(symbol, ensembl=None):
    q = urllib.parse.quote(symbol)
    links = {
        "BioLegend (Flow-Listed)": f"https://www.biolegend.com/en-us/search-results?Applications=FC&Keywords={q}",
        "CiteAb": f"https://www.citeab.com/antibodies/search?q={q}",
        "Biocompare": f"https://www.biocompare.com/Search-Antibodies/?search={q}&said=0&vcmpv=true",
    }
    if ensembl:
        links["Protein Atlas: Choroid Plexus"] = f"https://www.proteinatlas.org/{ensembl}-{q}/tissue/choroid+plexus"
    return links


def antibody_verdict(con, row):
    """Where an antibody's immunogen lies: from a stored epitope if mapped, else from its description."""
    if row["epitope_id"]:
        e = con.execute("SELECT start_aa, end_aa FROM epitope WHERE epitope_id = ?", (row["epitope_id"],)).fetchone()
        if e and e["start_aa"] and e["end_aa"]:
            verdict, _ = ranking.check_immunogen(con, row["target_id"], e["start_aa"], e["end_aa"])
            if verdict != "unknown":
                return verdict, f"residues {e['start_aa']}-{e['end_aa']} compared with the protein topology"
    return ranking.immunogen_from_text(con, row["target_id"], row["immunogen"])


def outside_segments(con, target_id):
    rows = con.execute("SELECT start_aa, end_aa FROM target_region WHERE target_id = ? AND region_type = 'extracellular' "
                       "AND start_aa IS NOT NULL ORDER BY start_aa", (target_id,)).fetchall()
    return ", ".join(f"{r['start_aa']}-{r['end_aa']}" for r in rows)


def log_result(con, target_id, symbol, choice, new, outcome, notes, ctx, species):
    """Save an antibody (if new) and one live-cell test result. Returns a message for the user."""
    warn = ""
    if choice is None:
        verdict, region = ranking.check_immunogen(con, target_id, new["start"], new["end"])
        epitope_id = None
        if region:
            con.execute("INSERT INTO epitope (target_id, epitope_type, region_type, start_aa, end_aa, mapping_method) "
                        "VALUES (?,?,?,?,?,?)", (target_id, "linear", region, new["start"], new["end"],
                                                 "immunogen residues compared with UniProt topology"))
            epitope_id = con.execute("SELECT last_insert_rowid()").fetchone()[0]
        status = "vendor_verified" if (new["vendor"] and new["catalog"]) else "needs_curation"
        con.execute("INSERT INTO antibody (target_id, epitope_id, clone, vendor, catalog_no, curation_status, reactivity) "
                    "VALUES (?,?,?,?,?,?,?)", (target_id, epitope_id, new["clone"] or None, new["vendor"] or None,
                                               new["catalog"] or None, status, species))
        antibody_id = con.execute("SELECT last_insert_rowid()").fetchone()[0]
        if verdict in ("intracellular", "removed"):
            warn = (" Warning: the immunogen you entered is not on the cell surface, so a live-cell result is unexpected. "
                    "Check that the cells were not permeabilized.")
    else:
        antibody_id = choice
    result = {"Worked": "works", "Did Not Work": "fails", "Mixed": "mixed"}[outcome]
    con.execute("INSERT INTO validation (antibody_id, application, cell_context, species, context_matches_cpec, "
                "nonpermeabilized, outcome, note) VALUES (?,?,?,?,?,?,?,?)",
                (antibody_id, "live_flow", ctx, species, 1, 1, result, notes or None))
    con.execute("UPDATE antibody SET stage = ? WHERE antibody_id = ?",
                ({"works": "passed", "fails": "failed", "mixed": "testing"}[result], antibody_id))
    if result == "works":
        con.execute("UPDATE antibody SET curation_status = 'lab_verified' WHERE antibody_id = ?", (antibody_id,))
    con.commit()
    return f"Saved: {symbol} test recorded as '{outcome}'.{warn}"


def save_registry_hit(con, target_id, hit, species):
    """Store a Registry record as a shortlisted antibody. Returns True if it was new."""
    if hit["rrid"]:
        row = con.execute("SELECT antibody_id FROM antibody WHERE target_id = ? AND rrid = ?",
                          (target_id, hit["rrid"])).fetchone()
        if row:
            con.execute("UPDATE antibody SET shortlisted = 1 WHERE antibody_id = ?", (row["antibody_id"],))
            con.commit()
            return False
    clonality = hit["clonality"] if hit["clonality"] in ("monoclonal", "polyclonal", "recombinant") else "unknown"
    con.execute(
        "INSERT INTO antibody (target_id, clone, host, isotype, clonality, reactivity, vendor, catalog_no, rrid, "
        "immunogen, conjugate, flow_claimed, applications_text, datasheet_url, citations, source, shortlisted, "
        "curation_status) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,1,?)",
        (target_id, hit["clone"] or None, hit["host"] or None, hit["isotype"] or None, clonality, species,
         hit["vendor"] or None, hit["catalog"] or None, hit["rrid"] or None, hit["epitope"] or None,
         hit["conjugate"] or None, "yes" if hit["flow_listed"] else "no", hit["applications"] or None,
         hit["vendor_url"] or None, hit["citations"], "Antibody Registry",
         "vendor_verified" if hit["vendor"] and hit["catalog"] else "needs_curation"))
    con.commit()
    return True


def export_zip(con):
    names = [r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type IN ('table','view') AND name NOT LIKE 'sqlite_%'")]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for n in names:
            if n == "setting":
                continue
            zf.writestr(f"{n}.csv", pd.read_sql_query(f"SELECT * FROM {n}", con).to_csv(index=False))
    return buf.getvalue()


def flag_counts(flags):
    return sum(1 for f in flags if f[0] == "problem"), sum(1 for f in flags if f[0] == "check")


con = db.connect()
st.title("CPEC Antibody Database")
st.caption("Finds and tracks antibodies that bind the outside of live cells, for sorting stem-cell-derived "
           "choroid plexus epithelial cells (CPECs).")
if "flash" in st.session_state:
    st.success(st.session_state.pop("flash"))

(tab_home, tab_find, tab_markers, tab_order, tab_track, tab_panel, tab_exp, tab_over, tab_records,
 tab_help) = st.tabs(["Test Next", "Find Antibodies", "Markers", "Order List", "Antibody Tracker", "Panel Designer",
                      "Experiments", "Overview", "My Records", "Help"])

# ================= TEST NEXT =================
with tab_home:
    top = st.columns([2, 3])
    if top[0].button("Update From Public Databases", type="primary",
                     help="Downloads membrane topology and facts from UniProt, choroid plexus RNA from the Human "
                          "Protein Atlas, and related papers from Europe PMC."):
        bar = st.progress(0.0, text="Starting")
        msgs = autofill.refresh_all(con, progress=lambda f, t: bar.progress(f, text=t))
        bar.empty()
        flash(msgs[0] if len(msgs) == 1 and "Could not reach" in msgs[0] else f"Updated {len(msgs)} markers.")
    top[1].caption("Click once after adding markers, then every few weeks. Needs an internet connection.")

    ranked = ranking.rank_targets(con)
    done = [r for r in ranked if r["proven"]]

    st.subheader(f"All {len(ranked)} Markers, Ranked")
    st.caption("Ranked automatically from live-cell access, how well the marker separates CPECs, the evidence on file, "
               "and papers found. Hover over a column name for its meaning.")
    table(pd.DataFrame([{
        "Rank": i, "Marker": r["symbol"], "Role": r["role_label"], "Score": r["score"],
        "Estimated Live-Cell Access": r["access"], "Longest Outside Segment (aa)": r["longest"] or None,
        "Specificity": r["specificity"], "Choroid Plexus RNA (nTPM)": r["cp_rna"], "Papers Found": r["n_papers"],
        "Antibodies on File": len(r["antibodies"]), "Evidence on File": r["evidence"]["text"]}
        for i, r in enumerate(ranked, 1)]).assign(**{"Choroid Plexus RNA (nTPM)": lambda d: num(d["Choroid Plexus RNA (nTPM)"])}))

    options = [r["symbol"] for r in ranked]
    first_open = next((i for i, r in enumerate(ranked) if not r["proven"]), 0)
    sym = st.selectbox("Open a Marker to Test", options, index=first_open, key="home_marker")
    r = ranked[options.index(sym)]
    st.markdown(f"### {r['symbol']}")
    st.caption(f"{L.tc(r['protein_name']) if r['protein_name'] else 'Name Not Loaded'} | {r['role_label']} | {r['klass']}")
    m = st.columns(5)
    m[0].metric("Estimated Live-Cell Access", r["access"])
    m[1].metric("Longest Outside Segment", f"{r['longest']} aa" if r["longest"] else "Not Loaded")
    m[2].metric("N-Glycosylation Sites", "Not Loaded" if r["n_glyc"] is None else r["n_glyc"])
    m[3].metric("Choroid Plexus RNA", f"{r['cp_rna']} nTPM" if r["cp_rna"] else "Not Loaded")
    m[4].metric("Papers Found", r["n_papers"])
    st.markdown(f"**Evidence on File:** {r['evidence']['text']}")
    with st.expander("Score Breakdown"):
        bd = pd.DataFrame(r["breakdown"])
        bd.loc[len(bd)] = {"Component": "Total", "Points": round(bd["Points"].sum(), 1), "Basis": ""}
        table(bd)
        st.caption("Points favor markers that are reachable on live cells, selective for CPECs, and backed by sources.")
    st.markdown("**Key Facts**")
    st.markdown("\n".join(f"- {x}" for x in r["reasons"]))
    if r["implications"]:
        st.markdown("**What to Know Before Choosing an Antibody**")
        st.markdown("\n".join(f"- {x}" for x in r["implications"]))
    st.markdown("  ".join(f"[{k}]({v})" for k, v in r["links"].items()))
    st.caption("Need an antibody? Open the Find Antibodies tab and search this marker.")
    if r["papers"]:
        st.markdown("**Papers to Check** (automatic leads, not verified)")
        for p_ in r["papers"][:5]:
            link = f"https://doi.org/{p_['doi']}" if p_["doi"] else None
            title = f"[{p_['title']}]({link})" if link else p_["title"]
            st.markdown(f"- {title} ({p_['journal']}, {p_['year']})")
    plan = ranking.test_plan(r["symbol"], r["role_label"], (r["antibodies"] or [None])[0], CELLS,
                             r["shed"], r["polarity"], r["contaminant"])
    with st.expander("Staining Plan With Controls"):
        st.code(plan, language=None)
        st.download_button("Download Plan (Text)", plan, file_name=f"plan_{r['symbol']}.txt", key=f"dl_{r['target_id']}")
    with st.form(f"log_{r['target_id']}"):
        st.markdown("**Record a Live-Cell Test** (cells not permeabilized)")
        names = [a["label"] for a in r["antibodies"]]
        other = "A Different Antibody (Fill In Below)"
        pick = st.selectbox("Antibody", names + [other], key=f"pick_{r['target_id']}")
        c = st.columns(3)
        vendor = c[0].text_input("Vendor", key=f"v_{r['target_id']}")
        catalog = c[1].text_input("Catalog Number", key=f"c_{r['target_id']}")
        clone = c[2].text_input("Clone", key=f"k_{r['target_id']}")
        with st.expander("Optional: Immunogen Residues (Checked Against the Topology)"):
            d = st.columns(2)
            istart = d[0].number_input("From Residue", min_value=0, value=0, key=f"s_{r['target_id']}")
            iend = d[1].number_input("To Residue", min_value=0, value=0, key=f"e_{r['target_id']}")
        outcome = st.radio("Result", ["Worked", "Did Not Work", "Mixed"], horizontal=True, index=None,
                           key=f"o_{r['target_id']}")
        cell_text = st.text_input("Cell Line and Differentiation Day (Optional, Shown in Your Test Record)",
                                  key=f"cl_{r['target_id']}")
        notes = st.text_input("Notes (Dilution, How Cells Were Detached)", key=f"n_{r['target_id']}")
        if st.form_submit_button("Save Result"):
            use_existing = pick != other
            new = {"vendor": vendor.strip(), "catalog": catalog.strip(), "clone": clone.strip(),
                   "start": int(istart), "end": int(iend)}
            if outcome is None:
                st.error("Choose a result.")
            elif not use_existing and not (new["vendor"] or new["clone"]):
                st.error("Enter at least a vendor or a clone for the new antibody.")
            else:
                choice = r["antibodies"][names.index(pick)]["id"] if use_existing else None
                ctx = CELLS + (f", {cell_text.strip()}" if cell_text.strip() else "")
                flash(log_result(con, r["target_id"], r["symbol"], choice, new, outcome, notes.strip(), ctx,
                                 r["species"]))

    st.subheader("Suggested Sort Strategy")
    table(pd.DataFrame(ranking.suggest_sort_strategy(ranked)))
    st.caption("The two best markers to keep and the two best leftover-stem-cell markers to remove, taken from the "
               "ranking above. The evidence column shows what is recorded for each marker, not a verdict.")
    if done:
        st.subheader("Worked on Live CPECs")
        table(pd.DataFrame([{"Marker": r_["symbol"], "Role": r_["role_label"],
                             "Antibody": r_["antibodies"][0]["label"] if r_["antibodies"] else ""} for r_ in done]))

# ================= FIND ANTIBODIES =================
with tab_find:
    st.write("Search the Antibody Registry for products against a marker. Each result shows vendor, catalog number, "
             "RRID, clone, conjugate and listed applications, with a datasheet check.")
    tg = con.execute("SELECT target_id, symbol, species, live_cell_accessible FROM target ORDER BY symbol").fetchall()
    c = st.columns([2, 1])
    sym = c[0].selectbox("Marker", [t["symbol"] for t in tg], key="find_marker")
    trow = [t for t in tg if t["symbol"] == sym][0]
    info = con.execute("SELECT * FROM target_info WHERE target_id = ?", (trow["target_id"],)).fetchone()
    segs = outside_segments(con, trow["target_id"])
    s_ = st.columns(3)
    s_[0].metric("Estimated Live-Cell Access", L.pretty(L.ACCESS, trow["live_cell_accessible"]))
    s_[1].metric("Outside Segments (Residues)", segs or "Not Loaded")
    s_[2].metric("N-Glycosylation Sites", "Not Loaded" if not info or info["n_glycosylation"] is None else info["n_glycosylation"])
    st.caption("Choose antibodies raised against residues inside the outside segments. An immunogen anywhere else "
               "cannot bind intact live cells.")
    if c[1].button("Search the Antibody Registry"):
        try:
            hits = autofill.registry_search(sym, trow["species"])
            st.session_state["reg"] = {"symbol": sym, "hits": hits}
        except autofill.SourceUnavailable:
            st.error("Could not reach the Antibody Registry. Check your internet connection and try again.")
    st.markdown("Also search: " + "  |  ".join(
        f"[{k}]({v})" for k, v in vendor_links(sym, info["hpa_ensembl"] if info else None).items()))

    reg = st.session_state.get("reg")
    if reg and reg["symbol"] == sym:
        for h in reg["hits"]:
            verdict, basis = ranking.immunogen_from_text(con, trow["target_id"], f"{h['epitope']} {h['name']}")
            h["verdict"], h["basis"] = verdict, basis
            h["flags"] = ranking.antibody_flags({**h, "azide_free": "unknown"}, verdict, basis)
            h["points"] = ranking.antibody_points({**h, "azide_free": "unknown"}, verdict)
        hits = sorted(reg["hits"], key=lambda x: -x["points"])
        st.write(f"{len(hits)} registry records for {sym}.")
        f = st.columns(3)
        only_flow = f[0].checkbox("Listed for Flow Cytometry Only")
        only_clonal = f[1].checkbox("Monoclonal or Recombinant Only")
        hide_inside = f[2].checkbox("Hide Immunogens Inside the Cell", value=True)
        if only_flow:
            hits = [h for h in hits if h["flow_listed"]]
        if only_clonal:
            hits = [h for h in hits if h["clonality"] in ("monoclonal", "recombinant")]
        if hide_inside:
            hits = [h for h in hits if h["verdict"] not in ("intracellular", "removed")]
        if hits:
            table(pd.DataFrame([{
                "Score": h["points"], "Vendor": h["vendor"], "Catalog": h["catalog"], "RRID": h["rrid"],
                "Name": h["name"], "Clone": h["clone"], "Clonality": L.pretty(L.CLONALITY, h["clonality"]),
                "Conjugate": h["conjugate"], "Isotype": h["isotype"], "Host": h["host"],
                "Listed Applications": h["applications"], "Citations": h["citations"],
                "Immunogen Location": L.pretty(L.VERDICT, h["verdict"]),
                "Problems": flag_counts(h["flags"])[0], "To Check": flag_counts(h["flags"])[1]} for h in hits]))
            st.caption("Score is a rough ordering: immunogen outside the cell, flow listing, conjugate, clonality, "
                       "citations. It does not measure live-cell performance.")
            options_ = {f"{h['vendor']} {h['catalog']} {h['name']}"[:110]: h for h in hits}
            look = st.selectbox("Datasheet Check for One Product", list(options_))
            for lvl, label, detail in options_[look]["flags"]:
                st.markdown(f"- **{TAG[lvl]} {label}.** {detail}")
            chosen = st.multiselect("Add Products to My Order List", list(options_))
            if st.button("Add to Order List") and chosen:
                n_new = sum(save_registry_hit(con, trow["target_id"], options_[k], trow["species"]) for k in chosen)
                flash(f"Added {len(chosen)} product(s) to the order list ({n_new} new).")
        else:
            st.info("No records match these filters.")

    with st.expander("Import Products From a CSV File"):
        st.caption("Columns it reads: Marker, Vendor, Catalog, Clone, Host, Isotype, Conjugate, Clonality, RRID. "
                   "Marker must already be in the database.")
        up = st.file_uploader("CSV File", type="csv", key="csv_up")
        if up is not None and st.button("Import Products"):
            df_ = pd.read_csv(up).fillna("")
            df_.columns = [c_.strip().lower() for c_ in df_.columns]
            n_ok, n_miss = 0, []
            for _, r_ in df_.iterrows():
                trg = con.execute("SELECT target_id, species FROM target WHERE UPPER(symbol) = ?",
                                  (str(r_.get("marker", "")).strip().upper(),)).fetchone()
                if not trg or not (str(r_.get("vendor", "")).strip() or str(r_.get("clone", "")).strip()):
                    n_miss.append(str(r_.get("marker", "?")))
                    continue
                cl = str(r_.get("clonality", "")).lower()
                con.execute("INSERT INTO antibody (target_id, clone, host, isotype, clonality, reactivity, vendor, catalog_no, rrid, "
                            "conjugate, source, shortlisted, curation_status) VALUES (?,?,?,?,?,?,?,?,?,?,?,1,?)",
                            (trg["target_id"], str(r_.get("clone", "")) or None, str(r_.get("host", "")) or None,
                             str(r_.get("isotype", "")) or None, cl if cl in ("monoclonal", "polyclonal", "recombinant") else "unknown",
                             trg["species"], str(r_.get("vendor", "")) or None, str(r_.get("catalog", "")) or None,
                             str(r_.get("rrid", "")) or None, str(r_.get("conjugate", "")) or None, "CSV import",
                             "vendor_verified" if str(r_.get("vendor", "")) and str(r_.get("catalog", "")) else "needs_curation"))
                n_ok += 1
            con.commit()
            flash(f"Imported {n_ok} product(s)." + (f" Skipped: {', '.join(n_miss)}." if n_miss else ""))

    with st.expander("Add a Product You Found Somewhere Else"):
        with st.form("manual_ab"):
            c1 = st.columns(3)
            v = c1[0].text_input("Vendor")
            cat = c1[1].text_input("Catalog Number")
            clone = c1[2].text_input("Clone")
            c2 = st.columns(4)
            host = c2[0].text_input("Host")
            iso = c2[1].text_input("Isotype")
            conj = c2[2].text_input("Conjugate (Blank if None)")
            clon = c2[3].selectbox("Clonality", ["Monoclonal", "Recombinant", "Polyclonal", "Unknown"])
            c3 = st.columns(2)
            azide = c3[0].selectbox("Azide-Free?", ["Unknown", "Yes", "No"])
            flow = c3[1].selectbox("Vendor Lists Flow Cytometry?", ["Unknown", "Yes", "No"])
            immuno = st.text_input("Immunogen From the Datasheet (For Example: Recombinant Protein, Residues 25-233)")
            url = st.text_input("Datasheet Link")
            if st.form_submit_button("Save Product"):
                if not (v.strip() or clone.strip()):
                    st.error("Enter at least a vendor or a clone.")
                else:
                    con.execute(
                        "INSERT INTO antibody (target_id, clone, host, isotype, clonality, reactivity, vendor, catalog_no, "
                        "immunogen, conjugate, azide_free, flow_claimed, datasheet_url, source, shortlisted, curation_status) "
                        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,1,?)",
                        (trow["target_id"], clone.strip() or None, host.strip() or None, iso.strip() or None,
                         clon.lower(), trow["species"], v.strip() or None, cat.strip() or None, immuno.strip() or None,
                         conj.strip() or None, azide.lower(), flow.lower(), url.strip() or None, "entered by hand",
                         "vendor_verified" if v.strip() and cat.strip() else "needs_curation"))
                    con.commit()
                    flash("Saved and added to the order list.")

# ================= MARKERS =================
with tab_markers:
    st.subheader("Discover More Markers")
    st.caption("Lists plasma membrane proteins with the highest choroid plexus RNA in the Human Protein Atlas. "
               "RNA is a lead, not proof of surface protein or of selectivity: many are also made in other tissues. "
               "Adding one downloads its topology, so you can see whether any part sits outside the cell.")
    if st.button("Find Candidates"):
        try:
            known = [t["symbol"] for t in con.execute("SELECT symbol FROM target").fetchall()]
            st.session_state["discover"] = autofill.discover_cp_markers(known)
        except autofill.SourceUnavailable:
            st.error("Could not reach the Human Protein Atlas. Check your internet connection and try again.")
    cands = st.session_state.get("discover")
    if cands is not None:
        if not cands:
            st.info("No new candidates found.")
        else:
            table(pd.DataFrame([{"Marker": c_["symbol"], "Choroid Plexus RNA (nTPM)": c_["cp_ntpm"],
                                 "Location (Human Protein Atlas)": c_["location"]} for c_ in cands]))
            chosen = st.multiselect("Add Markers", [c_["symbol"] for c_ in cands])
            if st.button("Add Selected Markers") and chosen:
                bar = st.progress(0.0, text="Adding")
                added, skipped = [], []
                for i, name in enumerate(chosen):
                    bar.progress(i / len(chosen), text=f"Adding {name}")
                    level = [c_ for c_ in cands if c_["symbol"] == name][0]["cp_ntpm"]
                    ok, msg = autofill.add_or_refresh_marker(
                        con, name, "human", "CPEC_identity",
                        note=f"Choroid plexus RNA {level} nTPM in the Human Protein Atlas. RNA lead only.")
                    (added if ok else skipped).append(name if ok else f"{name} ({msg})")
                bar.empty()
                st.session_state["discover"] = [c_ for c_ in cands if c_["symbol"] not in added]
                flash(f"Added {len(added)} marker(s)." + (f" Not added: {'; '.join(skipped)}" if skipped else ""))

    st.subheader("Add a Marker")
    st.caption("Type a gene symbol. The topology, facts and papers are downloaded and filled in automatically.")
    with st.form("add_marker"):
        c = st.columns([2, 1, 3])
        sym_new = c[0].text_input("Gene Symbol (For Example FOLR1)")
        sp_new = c[1].selectbox("Species", ["Human", "Mouse"])
        role = c[2].selectbox("What Is It For?", ROLE_OPTIONS)
        if st.form_submit_button("Add and Auto-Fill"):
            ok, msg = autofill.add_or_refresh_marker(con, sym_new, sp_new.lower(), ROLE_CODE[role])
            if ok:
                try:
                    row = con.execute("SELECT target_id, protein_name FROM target WHERE symbol = ? AND species = ?",
                                      (sym_new.strip(), sp_new.lower())).fetchone()
                    n = autofill.search_literature(con, row["target_id"], sym_new.strip(), row["protein_name"])
                    msg += f" {n} papers found."
                    if sp_new == "Human" and autofill.refresh_hpa(con, row["target_id"], sym_new.strip()):
                        msg += " Protein Atlas data added."
                except autofill.SourceUnavailable:
                    msg += " Could not reach Europe PMC or the Protein Atlas."
                flash(msg)
            else:
                st.error(msg)

    st.subheader("All Markers")
    rows = []
    for t in con.execute("SELECT t.*, i.n_glycosylation, i.hpa_cp_ntpm FROM target t "
                         "LEFT JOIN target_info i USING(target_id) ORDER BY t.symbol"):
        longest = con.execute("SELECT MAX(end_aa - start_aa + 1) FROM target_region WHERE target_id = ? "
                              "AND region_type = 'extracellular'", (t["target_id"],)).fetchone()[0]
        rows.append({"Marker": t["symbol"], "Protein": L.tc(t["protein_name"]) if t["protein_name"] else None, "Species": t["species"].capitalize(),
                     "Class": L.pretty(L.TARGET_CLASS, t["target_class"]),
                     "Estimated Live-Cell Access": L.pretty(L.ACCESS, t["live_cell_accessible"]),
                     "Longest Outside Segment (aa)": longest, "Specificity": L.pretty(L.SPECIFICITY, t["specificity"]),
                     "N-Glycosylation Sites": t["n_glycosylation"], "Choroid Plexus RNA (nTPM)": t["hpa_cp_ntpm"],
                     "UniProt Accession": t["uniprot_acc"], "UniProt Verified": "Yes" if t["uniprot_verified"] else "No"})
    tbl = pd.DataFrame(rows).assign(**{"Choroid Plexus RNA (nTPM)": lambda d: num(d["Choroid Plexus RNA (nTPM)"])})
    table(tbl)
    pick = st.selectbox("Open a Marker", tbl["Marker"].tolist(), key="open_marker")
    if pick:
        t = con.execute("SELECT * FROM target WHERE symbol = ?", (pick,)).fetchone()
        tid = t["target_id"]
        inf = con.execute("SELECT * FROM target_info WHERE target_id = ?", (tid,)).fetchone()
        longest = con.execute("SELECT MAX(end_aa - start_aa + 1) FROM target_region WHERE target_id = ? "
                              "AND region_type = 'extracellular'", (tid,)).fetchone()[0]
        pol = con.execute("SELECT polarity FROM cpec_context WHERE target_id = ? AND polarity != 'unknown'",
                          (tid,)).fetchone()
        st.markdown(f"### {pick}: {L.tc(t['protein_name']) if t['protein_name'] else ''}")
        st.write(t["accessibility_note"] or "")
        if inf:
            facts = [("What It Does", inf["function_text"]), ("Where It Sits in the Cell", inf["location_text"]),
                     ("Where It Is Expressed", inf["tissue_text"]), ("Modifications", inf["ptm_text"]),
                     ("N-Glycosylation Sites", inf["n_glycosylation"]), ("Disulfide Bonds", inf["n_disulfide"]),
                     ("Isoforms", inf["n_isoforms"]), ("Other Names", inf["aliases"]),
                     ("Protein Length (aa)", inf["length_aa"]), ("Mass (kDa)", inf["mass_kda"]),
                     ("Membrane-Spanning Segments", inf["n_tm"]), ("Solved Structures (PDB)", inf["n_pdb"]), ("Choroid Plexus RNA (nTPM, Human Protein Atlas)", inf["hpa_cp_ntpm"]),
                     ("Main Location (Human Protein Atlas)", inf["hpa_location"])]
            st.table(pd.DataFrame([(a, "" if b is None else str(b)) for a, b in facts if b not in (None, "")],
                                  columns=["Fact", "Value"]).set_index("Fact"))
        else:
            st.info("No facts downloaded yet. Click Update From Public Databases on the Test Next tab.")
        impl = ranking.marker_implications(t, inf, longest, pol["polarity"] if pol else None)
        if impl:
            st.markdown("**What This Means for Antibody Choice**")
            st.markdown("\n".join(f"- {x}" for x in impl))
        st.markdown("**Role in CPEC Work**")
        ctx_rows = con.execute("SELECT context, expected_in_cpec, polarity, evidence_basis, note FROM cpec_context "
                               "WHERE target_id = ?", (tid,)).fetchall()
        table(pd.DataFrame([{"Role": L.pretty(L.ROLE, r["context"]), "Expected in CPECs": L.pretty(L.EXPECTED, r["expected_in_cpec"]),
                             "Polarity": L.pretty(L.POLARITY, r["polarity"]),
                             "Evidence Basis": L.pretty(L.EVIDENCE_BASIS, r["evidence_basis"]), "Note": r["note"]}
                            for r in ctx_rows]))
        st.markdown("**Outside, Inside, and Membrane Segments**")
        seg_rows = con.execute("SELECT region_type, start_aa, end_aa, source FROM target_region WHERE target_id = ? "
                               "ORDER BY start_aa", (tid,)).fetchall()
        table(pd.DataFrame([{"Segment": L.pretty(L.REGION, r["region_type"]), "Start": r["start_aa"], "End": r["end_aa"],
                             "Length (aa)": (r["end_aa"] - r["start_aa"] + 1) if r["start_aa"] and r["end_aa"] else None,
                             "Source": r["source"]} for r in seg_rows]))
        mine = con.execute("SELECT * FROM antibody WHERE target_id = ?", (tid,)).fetchall()
        if mine:
            st.markdown("**Antibodies on File, With Datasheet Checks**")
            for a in mine:
                verdict, basis = antibody_verdict(con, a)
                flags = ranking.antibody_flags({**dict(a), "flow_listed": a["flow_claimed"]}, verdict, basis)
                bad, chk = flag_counts(flags)
                name = " ".join(x for x in [a["clone"] or "No Clone Name", a["vendor"], a["catalog_no"]] if x)
                with st.expander(f"{name} | {bad} Problems | {chk} To Check"):
                    for lvl, label, detail in flags:
                        st.markdown(f"- **{TAG[lvl]} {label}.** {detail}")

# ================= ORDER LIST =================
with tab_order:
    st.write("Products you shortlisted. Check each datasheet before ordering.")
    rows = con.execute("SELECT a.*, t.symbol FROM antibody a JOIN target t USING(target_id) "
                       "WHERE a.shortlisted = 1 ORDER BY t.symbol").fetchall()
    if not rows:
        st.info("Nothing here yet. Use Find Antibodies, then Add to Order List.")
    else:
        out = []
        for a in rows:
            verdict, basis = antibody_verdict(con, a)
            flags = ranking.antibody_flags({**dict(a), "flow_listed": a["flow_claimed"]}, verdict, basis)
            bad, chk = flag_counts(flags)
            out.append({"Marker": a["symbol"], "Vendor": a["vendor"], "Catalog": a["catalog_no"], "RRID": a["rrid"],
                        "Clone": a["clone"], "Conjugate": a["conjugate"],
                        "Clonality": L.pretty(L.CLONALITY, a["clonality"]),
                        "Immunogen Location": L.pretty(L.VERDICT, verdict), "Problems": bad, "To Check": chk,
                        "Datasheet": a["datasheet_url"], "_id": a["antibody_id"]})
        df = pd.DataFrame(out)
        table(df.drop(columns="_id"))
        st.download_button("Download Order List (CSV)", df.drop(columns="_id").to_csv(index=False),
                           file_name="order_list.csv")
        labels_ = {f"{r['Marker']} {r['Vendor']} {r['Catalog']} (ID {r['_id']})": r["_id"] for r in out}
        rm = st.multiselect("Remove From the List", list(labels_))
        if st.button("Remove Selected") and rm:
            for label in rm:
                con.execute("UPDATE antibody SET shortlisted = 0 WHERE antibody_id = ?", (labels_[label],))
            con.commit()
            flash("Removed.")

# ================= ANTIBODY TRACKER =================
with tab_track:
    st.write("Follow each antibody from candidate to tested: order and receipt dates, lot, price, titration and "
             "staining index. Saving a result on Test Next also moves the stage.")
    pool = con.execute("SELECT a.*, t.symbol FROM antibody a JOIN target t USING(target_id) "
                       "WHERE a.shortlisted = 1 OR a.stage != 'candidate' ORDER BY t.symbol").fetchall()
    if not pool:
        st.info("Nothing to track yet. Add products to the order list first.")
    else:
        pipe = pd.Series([L.pretty(L.STAGE, a["stage"]) for a in pool]).value_counts()
        mc = st.columns(len(L.STAGE))
        for i, (code, nm) in enumerate(L.STAGE.items()):
            mc[i].metric(nm, int(pipe.get(nm, 0)))
        table(pd.DataFrame([{"ID": a["antibody_id"], "Marker": a["symbol"], "Vendor": a["vendor"], "Catalog": a["catalog_no"],
                             "Clone": a["clone"], "Stage": L.pretty(L.STAGE, a["stage"]), "Ordered": a["ordered_on"],
                             "Received": a["received_on"], "Lot": a["lot"], "Price": a["price"]} for a in pool]))
        names = {f"{a['symbol']} | {a['vendor'] or 'No Vendor'} {a['catalog_no'] or a['clone'] or ''} (ID {a['antibody_id']})": a
                 for a in pool}
        with st.form("track_form"):
            who = st.selectbox("Antibody", list(names))
            c1 = st.columns(5)
            stg = c1[0].selectbox("Stage", list(L.STAGE.values()))
            od = c1[1].text_input("Date Ordered (YYYY-MM-DD)")
            rd = c1[2].text_input("Date Received (YYYY-MM-DD)")
            lot = c1[3].text_input("Lot Number")
            pr = c1[4].text_input("Price")
            if st.form_submit_button("Save Tracking"):
                code = {v: k for k, v in L.STAGE.items()}[stg]
                con.execute("UPDATE antibody SET stage = ?, ordered_on = COALESCE(NULLIF(?, ''), ordered_on), "
                            "received_on = COALESCE(NULLIF(?, ''), received_on), lot = COALESCE(NULLIF(?, ''), lot), "
                            "price = COALESCE(NULLIF(?, ''), price) WHERE antibody_id = ?",
                            (code, od.strip(), rd.strip(), lot.strip(), pr.strip(), names[who]["antibody_id"]))
                con.commit()
                flash("Tracking saved.")
        st.markdown("**Titration Log**")
        st.caption("Staining Index = (Positive Median - Negative Median) / (2 x Negative Standard Deviation). "
                   "The best dilution gives the highest index, not the brightest signal.")
        with st.form("titr_form"):
            who2 = st.selectbox("Antibody Titrated", list(names))
            t1 = st.columns(5)
            tdate = t1[0].text_input("Date (YYYY-MM-DD)")
            dil = t1[1].text_input("Dilution (For Example 1:100)")
            mp = t1[2].number_input("Positive Median", 0.0, value=0.0)
            mn = t1[3].number_input("Negative Median", 0.0, value=0.0)
            sd = t1[4].number_input("Negative Standard Deviation", 0.0, value=0.0)
            tn = st.text_input("Titration Notes")
            if st.form_submit_button("Save Titration Point"):
                if not dil.strip() or sd <= 0:
                    st.error("Enter a dilution and a negative standard deviation above zero.")
                else:
                    con.execute("INSERT INTO titration (antibody_id, run_date, dilution, mfi_pos, mfi_neg, sd_neg, note) "
                                "VALUES (?,?,?,?,?,?,?)", (names[who2]["antibody_id"], tdate.strip() or None, dil.strip(),
                                                          mp, mn, sd, tn.strip() or None))
                    con.commit()
                    flash("Titration point saved.")
        tit = con.execute("SELECT t.*, a.vendor, a.catalog_no, a.clone, g.symbol FROM titration t JOIN antibody a USING(antibody_id) "
                          "JOIN target g USING(target_id) ORDER BY g.symbol, t.antibody_id, t.titration_id").fetchall()
        if tit:
            rows_ = [{"Marker": r["symbol"], "Antibody": " ".join(x for x in [r["clone"], r["vendor"], r["catalog_no"]] if x),
                      "Date": r["run_date"], "Dilution": r["dilution"], "Positive Median": r["mfi_pos"],
                      "Negative Median": r["mfi_neg"], "Staining Index": round((r["mfi_pos"] - r["mfi_neg"]) / (2 * r["sd_neg"]), 1)}
                     for r in tit]
            tdf = pd.DataFrame(rows_)
            best = tdf.loc[tdf.groupby(["Marker", "Antibody"])["Staining Index"].idxmax()]
            st.markdown("**Best Dilution for Each Antibody**")
            table(best[["Marker", "Antibody", "Dilution", "Staining Index"]])
            table(tdf)

# ================= PANEL DESIGNER =================
with tab_panel:
    st.write("Pick products from your order list for the Keep and Remove gates. The check flags channel clashes, "
             "unconjugated antibodies, azide and immunogens that cannot bind live cells.")
    pool = con.execute("SELECT a.*, t.symbol FROM antibody a JOIN target t USING(target_id) WHERE a.shortlisted = 1 "
                       "ORDER BY t.symbol").fetchall()
    if not pool:
        st.info("Your order list is empty. Add products on the Find Antibodies tab first.")
    else:
        lab = {f"{a['symbol']} | {a['vendor'] or 'No Vendor'} {a['catalog_no'] or ''} | {a['conjugate'] or 'Unconjugated'}": a
               for a in pool}
        keep_sel = st.multiselect("Keep Gate (Markers on the Cells You Want)", list(lab))
        rem_sel = st.multiselect("Remove Gate (Markers on Cells You Want to Discard)", list(lab))
        members = []
        for role, sel in (("keep", keep_sel), ("remove", rem_sel)):
            for k in sel:
                a = lab[k]
                verdict, _ = antibody_verdict(con, a)
                members.append({"symbol": a["symbol"], "role": role, "conjugate": a["conjugate"],
                                "azide_free": a["azide_free"], "verdict": verdict, "label": k})
        if members:
            prow, pwarn = panels.check_panel(members)
            table(pd.DataFrame(prow))
            st.markdown("**Checks**")
            st.markdown("\n".join(f"- {w}" for w in pwarn))
            st.download_button("Download Panel (CSV)", pd.DataFrame(prow).to_csv(index=False), file_name="panel.csv")

# ================= EXPERIMENTS =================
with tab_exp:
    st.write("Log each sort or enrichment run. The tool totals purity gain, yield and viability across runs.")
    syms = [r[0] for r in con.execute("SELECT symbol FROM target ORDER BY symbol")]
    with st.form("exp_form"):
        e1 = st.columns(3)
        run_date = e1[0].text_input("Run Date (YYYY-MM-DD)")
        method = e1[1].selectbox("Method", ["FACS", "MACS", "Other"])
        cells = e1[2].text_input("Cell Line and Differentiation Day")
        keep_m = st.multiselect("Keep Markers Used", syms)
        rem_m = st.multiselect("Remove Markers Used", syms)
        e2 = st.columns(4)
        pb = e2[0].number_input("Purity Before (%)", 0.0, 100.0, 0.0)
        pa = e2[1].number_input("Purity After (%)", 0.0, 100.0, 0.0)
        yp = e2[2].number_input("Yield (%)", 0.0, 100.0, 0.0)
        vp = e2[3].number_input("Viability After (%)", 0.0, 100.0, 0.0)
        en = st.text_input("Notes")
        if st.form_submit_button("Save Run"):
            con.execute("INSERT INTO sort_experiment (run_date, keep_markers, remove_markers, cells, purity_before, "
                        "purity_after, yield_pct, viability_pct, method, note) VALUES (?,?,?,?,?,?,?,?,?,?)",
                        (run_date.strip() or None, ", ".join(keep_m), ", ".join(rem_m), cells.strip() or None,
                         pb, pa, yp, vp, method, en.strip() or None))
            con.commit()
            flash("Run saved.")
    runs = con.execute("SELECT * FROM sort_experiment ORDER BY run_date DESC, experiment_id DESC").fetchall()
    if runs:
        gains = [r["purity_after"] - r["purity_before"] for r in runs]
        m = st.columns(4)
        m[0].metric("Runs Logged", len(runs))
        m[1].metric("Average Purity Gain (Points)", round(sum(gains) / len(gains), 1))
        m[2].metric("Average Yield (%)", round(sum(r["yield_pct"] for r in runs) / len(runs), 1))
        m[3].metric("Average Viability (%)", round(sum(r["viability_pct"] for r in runs) / len(runs), 1))
        table(pd.DataFrame([{"Date": r["run_date"], "Method": r["method"], "Cells": r["cells"],
                             "Keep Markers": r["keep_markers"], "Remove Markers": r["remove_markers"],
                             "Purity Before (%)": r["purity_before"], "Purity After (%)": r["purity_after"],
                             "Yield (%)": r["yield_pct"], "Viability After (%)": r["viability_pct"],
                             "Notes": r["note"]} for r in runs]))
    else:
        st.info("No runs logged yet.")

# ================= OVERVIEW =================
with tab_over:
    st.write("What the database holds, what is missing, and what to do next.")
    q = lambda sql: con.execute(sql).fetchone()[0]
    o = st.columns(6)
    o[0].metric("Markers", q("SELECT COUNT(*) FROM target"))
    o[1].metric("Topology Loaded", q("SELECT COUNT(DISTINCT target_id) FROM target_region"))
    o[2].metric("Protein Atlas Data", q("SELECT COUNT(*) FROM target_info WHERE hpa_cp_ntpm IS NOT NULL"))
    o[3].metric("Antibody Records", q("SELECT COUNT(*) FROM antibody"))
    o[4].metric("On Order List", q("SELECT COUNT(*) FROM antibody WHERE shortlisted = 1"))
    o[5].metric("Your Tests", q("SELECT COUNT(*) FROM validation WHERE ref_id IS NULL"))
    st.markdown("**Markers by Role**")
    table(pd.DataFrame([{"Role": L.pretty(L.ROLE, r["context"]), "Markers": r["n"],
                         "With Topology": r["topo"], "With Antibody Records": r["abs"]} for r in con.execute(
        "SELECT c.context, COUNT(DISTINCT c.target_id) n, "
        "COUNT(DISTINCT CASE WHEN EXISTS (SELECT 1 FROM target_region g WHERE g.target_id = c.target_id) THEN c.target_id END) topo, "
        "COUNT(DISTINCT CASE WHEN EXISTS (SELECT 1 FROM antibody a WHERE a.target_id = c.target_id) THEN c.target_id END) abs "
        "FROM cpec_context c GROUP BY c.context")]))
    st.markdown("**Registry Coverage**")
    if st.button("Search the Registry for Every Marker"):
        humans = con.execute("SELECT target_id, symbol, species FROM target WHERE species = 'human' "
                             "AND target_class NOT IN ('glycan_epitope','glycoprotein_epitope')").fetchall()
        bar = st.progress(0.0, text="Searching")
        failed = 0
        for i, t in enumerate(humans):
            bar.progress(i / max(len(humans), 1), text=f"Searching {t['symbol']}")
            try:
                hits = autofill.registry_search(t["symbol"], t["species"])
            except autofill.SourceUnavailable:
                failed += 1
                if failed >= 3:
                    break
                continue
            vend = pd.Series([h["vendor"] for h in hits if h["vendor"]])
            con.execute("INSERT OR REPLACE INTO registry_summary VALUES (?,?,?,?,?,?,date('now'))",
                        (t["target_id"], len(hits), sum(h["flow_listed"] for h in hits),
                         sum(h["clonality"] in ("monoclonal", "recombinant") for h in hits),
                         sum(bool(h["conjugate"]) for h in hits), vend.mode().iloc[0] if len(vend) else None))
        con.commit()
        bar.empty()
        flash("Could not reach the Antibody Registry." if failed >= 3 else "Registry coverage updated.")
    cov = con.execute("SELECT t.symbol, s.* FROM registry_summary s JOIN target t USING(target_id) "
                      "ORDER BY s.n_flow DESC, t.symbol").fetchall()
    if cov:
        table(pd.DataFrame([{"Marker": r["symbol"], "Registry Records": r["n_records"], "Flow-Listed": r["n_flow"],
                             "Monoclonal or Recombinant": r["n_clonal"], "Conjugated": r["n_conjugated"],
                             "Most Common Vendor": r["top_vendor"], "Checked On": r["fetched_on"]} for r in cov]))
    else:
        st.caption("Not run yet. This counts vendor products per marker so you can see where antibodies are plentiful.")
    st.markdown("**Gaps to Fill**")
    gaps = []
    for t in con.execute("SELECT t.target_id, t.symbol, t.uniprot_acc, (SELECT COUNT(*) FROM target_region g WHERE g.target_id = t.target_id) topo, "
                         "(SELECT COUNT(*) FROM antibody a WHERE a.target_id = t.target_id AND a.vendor IS NOT NULL) vend, "
                         "(SELECT COUNT(*) FROM literature_hit h WHERE h.target_id = t.target_id) papers, "
                         "(SELECT hpa_cp_ntpm FROM target_info i WHERE i.target_id = t.target_id) hpa FROM target t ORDER BY t.symbol"):
        todo = []
        if not t["topo"]:
            todo.append("Load topology (Update From Public Databases)")
        if not t["vend"]:
            todo.append("Find a vendor product (Find Antibodies)")
        if not t["papers"]:
            todo.append("Search papers (Update From Public Databases)")
        if not t["hpa"]:
            todo.append("Load choroid plexus RNA")
        if todo:
            gaps.append({"Marker": t["symbol"], "Missing": len(todo), "Next Steps": "; ".join(todo)})
    table(pd.DataFrame(gaps).sort_values("Missing", ascending=False) if gaps else pd.DataFrame({"Marker": []}))

# ================= MY RECORDS =================
with tab_records:
    ab_rows = con.execute("SELECT a.*, t.symbol FROM antibody a JOIN target t USING(target_id) "
                          "ORDER BY t.symbol, a.clone").fetchall()
    recs = []
    for a in ab_rows:
        verdict, _ = antibody_verdict(con, a)
        live_apps = "('live_flow','facs_sort','macs','live_icc')"
        passes = con.execute(f"SELECT COUNT(*) FROM validation WHERE antibody_id = ? AND outcome = 'works' "
                             f"AND ref_id IS NULL AND application IN {live_apps}", (a["antibody_id"],)).fetchone()[0]
        reported = con.execute(f"SELECT COUNT(*) FROM validation WHERE antibody_id = ? AND outcome = 'works' "
                               f"AND ref_id IS NOT NULL AND application IN {live_apps}", (a["antibody_id"],)).fetchone()[0]
        recs.append({"ID": a["antibody_id"], "Marker": a["symbol"], "Clone": a["clone"], "Vendor": a["vendor"],
                     "Catalog": a["catalog_no"], "RRID": a["rrid"], "Host": a["host"], "Isotype": a["isotype"],
                     "Conjugate": a["conjugate"], "Clonality": L.pretty(L.CLONALITY, a["clonality"]),
                     "Immunogen Location": L.pretty(L.VERDICT, verdict),
                     "Flow Listed": L.pretty(L.YESNO, a["flow_claimed"]), "Azide-Free": L.pretty(L.YESNO, a["azide_free"]),
                     "Citations": a["citations"], "Status": L.pretty(L.STATUS, a["curation_status"]),
                     "Your Live-Cell Passes": passes, "Reported in Literature": reported, "Source": a["source"]})
    ab = pd.DataFrame(recs)
    f = st.columns([3, 2])
    text = f[0].text_input("Search Marker, Clone, or Vendor")
    hide = f[1].checkbox("Hide Unconfirmed (Clone Name Only)")
    if text and len(ab):
        mask = ab[["Marker", "Clone", "Vendor"]].fillna("").apply(lambda c: c.str.contains(text, case=False)).any(axis=1)
        ab = ab[mask]
    if hide and len(ab):
        ab = ab[ab["Status"] != "Clone Name Only"]
    table(ab)
    st.caption("Status: Clone Name Only means nothing else is confirmed. Details Missing means some fields are blank. "
               "Vendor Listed means vendor and catalog are recorded. Passed in Your Lab means a live-cell test worked.")
    tests = con.execute("SELECT t.symbol, a.clone, a.vendor, v.application, v.cell_context, v.outcome, v.note, "
                        "r.short_name FROM validation v JOIN antibody a USING(antibody_id) JOIN target t USING(target_id) "
                        "LEFT JOIN reference r ON r.ref_id = v.ref_id ORDER BY t.symbol").fetchall()

    def test_table(rs, source_col):
        return pd.DataFrame([{"Marker": r["symbol"], "Antibody": " ".join(x for x in [r["clone"], r["vendor"]] if x),
                              "Application": L.pretty(L.APPLICATION, r["application"]), "Cells": r["cell_context"],
                              "Result": L.pretty(L.OUTCOME, r["outcome"]), "Notes": r["note"],
                              source_col: r["short_name"] or "Your Test"} for r in rs])
    own_tests = [r for r in tests if r["short_name"] is None]
    lit_tests = [r for r in tests if r["short_name"] is not None]
    st.markdown("**Tests You Recorded**")
    if own_tests:
        table(test_table(own_tests, "Source"))
    else:
        st.info("You have not recorded any tests yet. Results you save on the Test Next tab appear here.")
    st.markdown("**Reported in the Literature (Not Your Tests)**")
    if lit_tests:
        table(test_table(lit_tests, "Reference"))
    else:
        st.caption("None on file.")
    with st.expander("Fix Vendor, Catalog Number, or RRID"):
        if len(ab):
            with st.form("fix_ab"):
                ids = ab["ID"].tolist()
                sel = st.selectbox("Antibody", ids, format_func=lambda i: f"{ab[ab['ID']==i]['Marker'].iloc[0]} "
                                                                           f"{ab[ab['ID']==i]['Clone'].iloc[0]} (ID {i})")
                c = st.columns(3)
                v, cat, rr = c[0].text_input("Vendor"), c[1].text_input("Catalog Number"), c[2].text_input("RRID")
                if st.form_submit_button("Save"):
                    sets, vals = [], []
                    for col, val in [("vendor", v), ("catalog_no", cat), ("rrid", rr)]:
                        if val.strip():
                            sets.append(f"{col} = ?")
                            vals.append(val.strip())
                    if v.strip() and cat.strip():
                        sets.append("curation_status = CASE WHEN curation_status IN ('seed_unverified','needs_curation') "
                                    "THEN 'vendor_verified' ELSE curation_status END")
                    if sets:
                        con.execute(f"UPDATE antibody SET {', '.join(sets)} WHERE antibody_id = ?", (*vals, int(sel)))
                        con.commit()
                        flash("Saved.")
                    else:
                        st.warning("Nothing to save.")
    st.markdown("**Download Your Data**")
    d = st.columns(2)
    d[0].download_button("All Tables as CSV (Zip)", export_zip(con), file_name="cpec_antibody_db_csv.zip")
    d[1].download_button("Database File", db.DB.read_bytes(), file_name="cpec_antibodies.sqlite")

# ================= HELP =================
with tab_help:
    st.markdown("""
**What This Does.** To sort living cells, an antibody must bind a part of the protein that sticks out of the cell.
This database finds out which parts stick out (UniProt), ranks which markers to test, finds products (Antibody
Registry), checks their datasheets, and keeps your test results.

**Your Routine**
1. **Test Next:** click *Update From Public Databases*, then open the top marker.
2. **Find Antibodies:** search the Registry for that marker, read the datasheet checks, add products to the order list.
3. Run the staining plan on live, unpermeabilized cells.
4. Back on **Test Next**, record *Worked*, *Did Not Work*, or *Mixed*.

**Glossary**
- **Estimated Live-Cell Access:** a rule of thumb from UniProt topology. Good means a GPI anchor, or an outside
  segment of 30 amino acids or more (single-pass proteins) or 60 or more (multi-pass proteins). Partial means only
  short outside loops. None means no outside segment.
- **Longest Outside Segment (aa):** the longest stretch of the protein outside the cell, in amino acids.
- **Immunogen:** the piece of protein used to raise the antibody. It must lie outside the cell.
- **Specificity:** Selective means mostly CPECs. Shared With Other Cells and Found Everywhere mean the marker needs a partner.
- **nTPM:** normalized transcripts per million, the RNA level reported by the Human Protein Atlas.
- **RRID:** a permanent identifier for an antibody product that you can cite in papers.
- **Conjugate:** a fluorophore attached to the antibody. Unconjugated antibodies need a labeled secondary.
- **Azide-Free:** buffer without sodium azide, which is toxic to cells you plan to culture after sorting.
- **Status:** Clone Name Only, Details Missing, Vendor Listed, or Passed in Your Lab.

**How the Score Works.** Each marker earns points for live-cell access (up to 3), selectivity (from minus 2 to 1),
its role in the sort (up to 2), evidence for the CPEC link (up to 1), and papers found (up to 2). Open the Score
Breakdown on any marker to see its own numbers.

**Reading the Datasheet Checks**
- *Immunogen Outside the Cell* is required. An antibody raised against an inside region cannot bind intact live cells.
- *Listed for Flow Cytometry* shows the vendor tested flow, not that they used live, unfixed cells.
- *Azide* is toxic to cells you plan to culture after sorting.
- Unconjugated antibodies need a matched secondary and a secondary-only control.

**Limits.** The Registry lists what vendors submit, so entries can be incomplete. Azide status is not in the
Registry, so it shows as unknown until you check the datasheet. Scores and access estimates are rules of
thumb; only your own test settles whether an antibody works on your cells. Product details change, so confirm
catalog numbers and prices with the vendor before ordering.
""")
con.close()
