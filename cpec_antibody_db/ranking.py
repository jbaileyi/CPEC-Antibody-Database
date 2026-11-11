"""Scoring and recommendations. Everything here is computed from the data, nothing is typed in by hand."""
import re
import urllib.parse

import labels

LIVE_APPS = ("live_flow", "facs_sort", "macs", "live_icc")
CURATION_RANK = {"lab_verified": 3, "vendor_verified": 2, "needs_curation": 1, "seed_unverified": 0}
ACCESS_SCORE = {"yes": 3.0, "partial": 1.5, "unknown": 0.5, "no": -5.0}
SPECIFICITY_SCORE = {"selective": 1.0, "shared": -0.5, "ubiquitous": -2.0, "unknown": 0.0}
ROLE_LABEL = labels.ROLE
ROLE_SCORE = {"CPEC_identity": 2.0, "pluripotent_contaminant": 1.0, "neural_progenitor_contaminant": 1.0,
              "endothelial_contaminant": 0.5, "mesenchymal_contaminant": 0.5, "pan_epithelial": 0.5}
EVIDENCE_SCORE = {"search_verified": 1.0, "prior_knowledge": 0.5, "hypothesis": 0.0}


def check_immunogen(con, target_id, start, end):
    """Compare an antibody's immunogen residues with the protein topology.
    Returns (verdict, region_type) where verdict is extracellular, intracellular, mixed, or unknown."""
    if not start or not end or end < start:
        return "unknown", None
    regions = con.execute("SELECT region_type, start_aa, end_aa, note FROM target_region WHERE target_id = ? "
                          "AND start_aa IS NOT NULL AND end_aa IS NOT NULL", (target_id,)).fetchall()
    if not regions:
        return "unknown", None
    total = end - start + 1
    ext = other = removed = 0
    for r in regions:
        overlap = max(0, min(end, r["end_aa"]) - max(start, r["start_aa"]) + 1)
        if r["region_type"] == "extracellular":
            ext += overlap
        elif r["region_type"] in ("cytoplasmic", "transmembrane", "lumenal"):
            other += overlap
        elif r["region_type"] == "other" and "removed" in (r["note"] or "").lower():
            removed += overlap
    if ext / total >= 0.8:
        return "extracellular", "extracellular"
    if ext == 0 and other == 0 and removed > 0:
        return "removed", "other"
    if ext == 0 and other > 0:
        return "intracellular", "cytoplasmic"
    if ext > 0:
        return "mixed", "other"
    return "unknown", None


def evidence_summary(con, target_id, role):
    """What is actually recorded for a marker. Nothing is assumed: the text names the kind of evidence on file."""
    rows = con.execute(
        "SELECT v.outcome, v.application, v.context_matches_cpec, v.cell_context, v.ref_id, r.short_name "
        "FROM validation v JOIN antibody a USING(antibody_id) LEFT JOIN reference r ON r.ref_id = v.ref_id "
        "WHERE a.target_id = ?", (target_id,)).fetchall()
    own = [r for r in rows if r["ref_id"] is None]
    if own:
        n = {k: sum(1 for r in own if r["outcome"] == k) for k in ("works", "fails", "mixed")}
        return {"level": "own", "text": f"Your tests: {n['works']} worked, {n['fails']} did not work, {n['mixed']} mixed"}
    reported = [r for r in rows if r["outcome"] == "works" and r["application"] in LIVE_APPS]
    if reported:
        r = reported[0]
        where = "CPEC-directed cells" if r["context_matches_cpec"] else (r["cell_context"] or "another cell type")
        return {"level": "reported", "text": f"Live-cell sort reported in {where} ({r['short_name']}); not tested in CPECs here"}
    basis = {"search_verified": "Supported by a retrieved source (not a live-cell test)",
             "prior_knowledge": "General knowledge only, no test or source recorded",
             "hypothesis": "Hypothesis only, no test or source recorded"}
    return {"level": role["evidence_basis"], "text": basis[role["evidence_basis"]]}


def _links(symbol, accession):
    out = {"PubMed": "https://pubmed.ncbi.nlm.nih.gov/?term=" + urllib.parse.quote(f'"choroid plexus" {symbol}'),
           "Protein Atlas": f"https://www.proteinatlas.org/search/{urllib.parse.quote(symbol)}"}
    if accession:
        out["UniProt"] = f"https://www.uniprot.org/uniprotkb/{accession}/entry"
    return out


def rank_targets(con):
    """Return one dict per target, best candidates first. Proven markers are flagged, not removed."""
    out = []
    for t in con.execute("SELECT * FROM target").fetchall():
        contexts = con.execute("SELECT context, expected_in_cpec, evidence_basis, polarity FROM cpec_context WHERE target_id = ?",
                               (t["target_id"],)).fetchall()
        if not contexts:
            continue
        role = max(contexts, key=lambda c: ROLE_SCORE.get(c["context"], 0))
        longest = con.execute("SELECT MAX(end_aa - start_aa + 1) FROM target_region "
                              "WHERE target_id = ? AND region_type = 'extracellular'", (t["target_id"],)).fetchone()[0]
        info = con.execute("SELECT * FROM target_info WHERE target_id = ?", (t["target_id"],)).fetchone()
        papers = con.execute("SELECT title, year, journal, doi, source_id FROM literature_hit WHERE target_id = ? "
                             "ORDER BY year DESC", (t["target_id"],)).fetchall()
        proven = con.execute(
            "SELECT 1 FROM validation v JOIN antibody a USING(antibody_id) WHERE a.target_id = ? AND v.outcome = 'works' "
            "AND v.application IN (?,?,?,?) AND v.context_matches_cpec = 1 LIMIT 1", (t["target_id"], *LIVE_APPS)).fetchone()
        abs_ = []
        n_intracellular = 0
        for a in con.execute(
                "SELECT a.antibody_id, a.clone, a.vendor, a.catalog_no, a.curation_status, a.isotype, a.conjugate, e.region_type "
                "FROM antibody a LEFT JOIN epitope e USING(epitope_id) WHERE a.target_id = ?", (t["target_id"],)):
            if a["region_type"] in ("cytoplasmic", "transmembrane"):
                n_intracellular += 1
                continue
            label = " ".join(x for x in [a["clone"] or "(no clone name)", a["vendor"], a["catalog_no"]] if x)
            abs_.append({"id": a["antibody_id"], "label": label, "status": a["curation_status"],
                         "isotype": a["isotype"], "conjugate": a["conjugate"]})
        abs_.sort(key=lambda x: (-CURATION_RANK.get(x["status"], 0), x["id"]))

        parts = [
            ("Live-Cell Access", ACCESS_SCORE.get(t["live_cell_accessible"], 0),
             f"Estimated access is {labels.pretty(labels.ACCESS, t['live_cell_accessible'])}"),
            ("Selectivity", SPECIFICITY_SCORE.get(t["specificity"], 0),
             labels.pretty(labels.SPECIFICITY, t["specificity"])),
            ("Role in the Sort", ROLE_SCORE.get(role["context"], 0), labels.ROLE.get(role["context"], role["context"])),
            ("Evidence for the CPEC Link", EVIDENCE_SCORE.get(role["evidence_basis"], 0),
             labels.pretty(labels.EVIDENCE_BASIS, role["evidence_basis"])),
            ("Papers Found", 0.5 * min(len(papers), 4), f"{len(papers)} found (up to 4 count)"),
        ]
        score = sum(pt for _, pt, _ in parts)
        reasons = []
        if longest:
            reasons.append(f"Outside the cell: longest exposed segment is {longest} aa "
                           f"(access estimate: {t['live_cell_accessible']}).")
        else:
            reasons.append(f"Access estimate: {t['live_cell_accessible']} (no topology loaded yet; click Update).")
        if t["specificity"] in ("shared", "ubiquitous"):
            reasons.append(f"Also found on other cell types ({t['specificity']}), so it will not isolate CPECs alone.")
        reasons.append({"search_verified": "Link to CPECs is supported by a retrieved source.",
                        "prior_knowledge": "Link to CPECs comes from general knowledge and is not re-verified.",
                        "hypothesis": "Link to CPECs is a hypothesis."}[role["evidence_basis"]])
        if papers:
            reasons.append(f"{len(papers)} paper(s) mention it with choroid plexus and sorting or surface terms.")
        if n_intracellular:
            reasons.append(f"{n_intracellular} antibody on file targets an inside-the-cell region and is hidden.")
        out.append({
            "target_id": t["target_id"], "symbol": t["symbol"], "protein_name": t["protein_name"],
            "species": t["species"], "evidence": evidence_summary(con, t["target_id"], role),
            "context": role["context"], "role_label": ROLE_LABEL.get(role["context"], role["context"]),
            "score": round(score, 2), "proven": bool(proven), "reasons": reasons,
            "breakdown": [{"Component": n, "Points": round(pt, 1), "Basis": b} for n, pt, b in parts],
            "access": labels.pretty(labels.ACCESS, t["live_cell_accessible"]), "access_code": t["live_cell_accessible"],
            "specificity": labels.pretty(labels.SPECIFICITY, t["specificity"]), "longest": longest or 0,
            "cp_rna": (info["hpa_cp_ntpm"] if info and info["hpa_cp_ntpm"] not in (None, "") else ""),
            "n_glyc": (info["n_glycosylation"] if info else None), "n_papers": len(papers),
            "klass": labels.pretty(labels.TARGET_CLASS, t["target_class"]),
            "best_status": labels.pretty(labels.STATUS, abs_[0]["status"]) if abs_ else "None on File",
            "antibodies": abs_, "papers": [dict(p) for p in papers],
            "links": _links(t["symbol"], t["uniprot_acc"]),
            "implications": marker_implications(t, info, longest, role["polarity"]),
            "shed": bool(info and info["shedding_hint"]), "polarity": role["polarity"],
            "contaminant": role["context"].endswith("_contaminant"),
            "action": "Test an antibody you have or buy" if abs_ else "Find an antibody (none on file)",
        })
    out.sort(key=lambda r: (-r["score"], r["symbol"]))
    return out


def suggest_sort_strategy(ranked):
    """Draft panel: two best markers to keep, two best leftover-stem-cell markers to remove."""
    keep = [r for r in ranked if r["context"] == "CPEC_identity"][:2]
    remove = [r for r in ranked if r["context"] == "pluripotent_contaminant"][:2]
    rows = []
    for r, sel in [(x, "keep") for x in keep] + [(x, "remove") for x in remove]:
        rows.append({"Marker": r["symbol"], "Action": sel.capitalize(),
                     "Antibody": r["antibodies"][0]["label"] if r["antibodies"] else "None on File Yet",
                     "Evidence": r["evidence"]["text"]})
    return rows


# ---------- Antibody datasheet checks ----------
def immunogen_from_text(con, target_id, text):
    """Judge where an antibody's immunogen lies, from a datasheet or registry description.
    Returns (verdict, basis). Residue ranges are checked against the topology; keywords are a weaker fallback."""
    text = text or ""
    m = re.search(r"(\d{1,4})\s*(?:-|\u2013|to)\s*(\d{1,4})", text)
    if m and int(m.group(1)) < int(m.group(2)):
        verdict, _ = check_immunogen(con, target_id, int(m.group(1)), int(m.group(2)))
        if verdict != "unknown":
            return verdict, f"residues {m.group(1)}-{m.group(2)} compared with the protein topology"
    low = text.lower()
    if re.search(r"extracellular|ectodomain|extra-cellular", low):
        return "extracellular", "described as extracellular"
    if re.search(r"cytoplasmic|intracellular|cytosolic", low):
        return "intracellular", "described as inside the cell"
    return "unknown", "not stated" if not text.strip() else "terminus named but location not clear; check the topology"


def antibody_flags(ab, verdict, basis):
    """Datasheet checklist for one antibody. Each item is (level, label, detail); level is ok, check, or problem."""
    f = []
    if verdict == "extracellular":
        f.append(("ok", "Immunogen Outside the Cell", basis))
    elif verdict == "intracellular":
        f.append(("problem", "Immunogen Inside the Cell", "It cannot bind intact live cells. " + basis))
    elif verdict == "removed":
        f.append(("problem", "Immunogen in a Removed Segment",
                  "That stretch is cut off when the GPI anchor is added, so it is not on the cell. " + basis))
    elif verdict == "mixed":
        f.append(("check", "Immunogen Spans Inside and Outside", "Only part of it is reachable; " + basis))
    else:
        f.append(("check", "Immunogen Location Unknown", "Read the datasheet for the immunogen and compare it with the outside segments."))
    flow = ab.get("flow_listed", ab.get("flow_claimed"))
    if flow in (True, "yes"):
        f.append(("ok", "Listed for Flow Cytometry", "A flow listing does not say the test used live, unfixed cells. Confirm on the datasheet."))
    elif flow in (False, "no"):
        f.append(("check", "Not Listed for Flow Cytometry", "You would be the first to validate it for flow."))
    else:
        f.append(("check", "Flow Listing Unknown", "Check the vendor datasheet."))
    az = ab.get("azide_free", "unknown")
    if az == "yes":
        f.append(("ok", "Azide-Free", "Suitable for sorting cells you will culture."))
    elif az == "no":
        f.append(("problem", "Contains Sodium Azide", "Azide is toxic to cells you plan to culture after sorting. Wash well or choose an azide-free version."))
    else:
        f.append(("check", "Azide Status Unknown", "For sorted cells you will culture, look for azide-free formats (often labeled LEAF or NA/LE)."))
    conj = (ab.get("conjugate") or "").strip().lower()
    if conj and not re.search(r"unconjugated|none|purified|^no$", conj):
        f.append(("ok", f"Conjugated ({ab.get('conjugate')})", "One-step staining."))
    else:
        f.append(("check", "Unconjugated or Conjugate Unknown", "Needs a labeled secondary matched to host and isotype, plus a secondary-only control."))
    cl = (ab.get("clonality") or "").lower()
    if cl in ("monoclonal", "recombinant"):
        f.append(("ok", cl.capitalize(), "Consistent from lot to lot."))
    elif cl == "polyclonal":
        f.append(("check", "Polyclonal", "Varies between lots and is seldom ideal for sorting."))
    n = ab.get("citations") or 0
    f.append(("ok" if n else "check", f"Cited in {n} Papers" if n else "No Citations Listed", "Citations show use, not live-cell success."))
    return f


def antibody_points(ab, verdict):
    """Rough ordering for candidate lists. A heuristic, not a measurement."""
    pts = {"extracellular": 3, "intracellular": -6, "removed": -6}.get(verdict, 0)
    pts += 2 if ab.get("flow_listed") in (True, "yes") else 0
    pts += {"yes": 1, "no": -1}.get(ab.get("azide_free", "unknown"), 0)
    conj = (ab.get("conjugate") or "").strip().lower()
    pts += 1 if conj and not re.search(r"unconjugated|none|purified|^no$", conj) else 0
    pts += 1 if (ab.get("clonality") or "").lower() in ("monoclonal", "recombinant") else 0
    pts += 0.4 * min(ab.get("citations") or 0, 5)
    return round(pts, 1)


def marker_implications(t, info, longest, polarity):
    """What the marker's facts mean for choosing and using an antibody."""
    out = []
    if t["target_class"] == "multi_pass_tm" and longest and longest < 60:
        out.append(f"Only short outside loops (longest {longest} aa). Many catalog antibodies for proteins like this "
                   "target the inside tail, so confirm the immunogen before buying.")
    if info:
        if (info["n_glycosylation"] or 0) >= 3:
            out.append(f"{info['n_glycosylation']} N-glycosylation sites. Sugars can hide or form epitopes, so an antibody may "
                       "bind differently in your cells than in the vendor's test cells.")
        if info["shedding_hint"]:
            out.append("UniProt text mentions shedding or cleavage. Soluble fragments in the medium can soak up antibody, "
                       "so wash cells before staining.")
        if (info["n_isoforms"] or 0) >= 2:
            out.append(f"{info['n_isoforms']} isoforms. Check that the immunogen is present in the isoform your cells express.")
        if (info["n_disulfide"] or 0) >= 2:
            out.append("Disulfide-stabilized folds mean conformational epitopes: avoid reducing agents in the staining buffer.")
        if info["hpa_cp_ntpm"] not in (None, ""):
            out.append(f"Choroid plexus RNA in the Human Protein Atlas: {info['hpa_cp_ntpm']} nTPM (RNA level, not surface protein).")
    if polarity in ("apical", "basolateral"):
        out.append(f"Reported {polarity} in choroid plexus epithelium. In intact cysts or monolayers with tight junctions "
                   "the antibody may not reach it, so stain dissociated cells.")
    if t["specificity"] in ("shared", "ubiquitous"):
        out.append(f"Found on other cell types ({t['specificity']}), so pair it with a second marker.")
    return out


def test_plan(symbol, role_label, ab, cell_context, shed, polarity, contaminant):
    """A staining plan with controls, in plain text. General practice, to adjust to the user's protocol."""
    ab = ab or {}
    iso = ab.get("isotype") or "matched"
    conj = ab.get("conjugate") or "the same conjugate"
    neg = ("your CPEC-directed cells (which should be negative or low)" if contaminant
           else "undifferentiated stem cells (which should be negative or low for most CPEC markers)")
    name = ab.get("label") or "the antibody"
    lines = [f"Live-Cell Staining Plan: {symbol} ({role_label})", f"Cells: {cell_context}", f"Antibody: {name}", "",
             "Controls (Run All)",
             "1. Unstained cells.",
             f"2. Isotype control ({iso}, {conj}) or a fluorescence-minus-one control.",
             "3. Viability dye (DAPI, 7-AAD or propidium iodide). Dead cells bind antibody nonspecifically; gate them out.",
             f"4. Negative reference: {neg}.",
             "5. Secondary-only control, if the primary antibody is unconjugated.", "",
             "Steps",
             "1. Detach gently with Accutase or EDTA rather than trypsin. Trypsin can cut extracellular domains and remove the epitope.",
             "2. If an enzyme was used, rest the cells 30 to 60 minutes in medium before staining.",
             "3. Filter to a single-cell suspension, count, and keep cells at 4 C from here on. Cold slows antibody-induced internalization.",
             "4. Wash in cold buffer without sodium azide (PBS with 0.5 to 2% BSA and 2 mM EDTA).",
             "5. Optional: block Fc receptors. Choroid plexus stroma contains macrophages that bind antibody Fc regions.",
             "6. Titrate the antibody (three dilutions) before the real experiment. Stain 20 to 30 minutes on ice, in the dark. Do NOT fix or permeabilize first.",
             "7. Wash twice, add the viability dye, and acquire promptly.",
             "8. Gate: cells, then single cells, then live cells, then the marker. Set the positive gate with the controls above.", "",
             "RECORD 'WORKED' only if a clear positive population separates from the isotype or FMO control and from the negative reference.",
             "Record 'Did not work' if the signal matches the controls, and 'Mixed' if separation is weak or varies between runs."]
    cautions = []
    if shed:
        cautions.append("This protein may be shed, so wash cells before staining.")
    if polarity in ("apical", "basolateral"):
        cautions.append(f"It is reported {polarity}; stain dissociated cells, not intact cysts.")
    if cautions:
        lines += ["", "Cautions for This Marker"] + [f"- {c}" for c in cautions]
    return "\n".join(lines)
