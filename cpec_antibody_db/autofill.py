"""Automatic data collection from public sources.

UniProt supplies membrane topology (which parts of a protein are outside the cell).
Europe PMC supplies papers that mention the marker together with choroid plexus and sorting terms.
Neither call needs an account or key. Nothing here talks to antibody vendors, because they offer no
open interface; vendor details are entered by the user.
"""
import json
import re
import time
import urllib.parse
import urllib.request
from datetime import date

UA = {"User-Agent": "cpec-antibody-db/0.2 (research database)"}
TAXID = {"human": 9606, "mouse": 10090}
UNIPROT_FIELDS = ("accession,gene_primary,protein_name,length,ft_topo_dom,ft_transmem,ft_signal,ft_lipid,"
                  "ft_carbohyd,ft_disulfid,cc_function,cc_subcellular_location,cc_tissue_specificity,"
                  "cc_ptm,cc_alternative_products,gene_names,mass,xref_pdb")
SURFACE_CLASSES = ("single_pass_tm", "multi_pass_tm", "gpi_anchored")


class SourceUnavailable(Exception):
    """Raised when a public database cannot be reached."""


class NotSurfaceProtein(Exception):
    """Raised when UniProt shows no membrane anchor, so the protein cannot be a live-cell surface target."""


def _get_json(url, timeout=25):
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.load(resp)
    except Exception as exc:
        raise SourceUnavailable(str(exc)) from exc


# ---------- UniProt ----------
def fetch_uniprot(symbol, species="human"):
    taxid = TAXID.get(species)
    if taxid is None:
        return None
    query = f"gene_exact:{symbol} AND organism_id:{taxid} AND reviewed:true"
    url = ("https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(query)
           + f"&fields={UNIPROT_FIELDS}&format=json&size=1")
    results = _get_json(url).get("results", [])
    return results[0] if results else None


def parse_entry(entry):
    """Turn a UniProt JSON entry into a plain dict."""
    genes = entry.get("genes") or []
    gene = genes[0].get("geneName", {}).get("value") if genes else None
    desc = entry.get("proteinDescription", {})
    name = (desc.get("recommendedName", {}).get("fullName", {}).get("value")
            or (desc.get("submissionNames") or [{}])[0].get("fullName", {}).get("value"))
    segments, n_tm, gpi, signal, n_glyc, n_ss, gpi_site = [], 0, False, False, 0, 0, None
    for f in entry.get("features", []):
        loc = f.get("location", {})
        start, end = loc.get("start", {}).get("value"), loc.get("end", {}).get("value")
        ftype, fdesc = f.get("type", ""), (f.get("description") or "")
        if ftype == "Transmembrane":
            n_tm += 1
            segments.append(("transmembrane", start, end, fdesc))
        elif ftype == "Topological domain":
            low = fdesc.lower()
            kind = ("extracellular" if low.startswith("extracellular") else
                    "cytoplasmic" if low.startswith("cytoplasmic") else
                    "lumenal" if low.startswith("lumenal") else "other")
            segments.append((kind, start, end, fdesc))
        elif ftype == "Signal":
            signal = True
        elif ftype == "Lipidation" and "GPI" in fdesc.upper():
            gpi, gpi_site = True, start
        elif ftype == "Glycosylation" and "N-linked" in fdesc:
            n_glyc += 1
        elif ftype == "Disulfide bond":
            n_ss += 1
    texts = _comment_texts(entry)
    isoforms = sum(len(c.get("isoforms", [])) for c in entry.get("comments", [])
                   if c.get("commentType") == "ALTERNATIVE PRODUCTS")
    shed = bool(re.search(r"\b(shed|shedding|ectodomain|soluble form|released from the cell|proteolytic(ally)? cleav)",
                          (texts["ptm"] + " " + texts["function"]).lower()))
    syn = [x.get("value") for g in genes for x in (g.get("synonyms") or []) if x.get("value")]
    alts = [a.get("fullName", {}).get("value") for a in desc.get("alternativeNames", []) or []]
    alts += [a.get("fullName", {}).get("value") for a in desc.get("cdAntigenNames", []) or []]
    alts += [a.get("value") for a in desc.get("cdAntigenNames", []) or [] if isinstance(a, dict) and a.get("value")]
    aliases = []
    for a in syn + alts:
        if a and a not in aliases and a != name:
            aliases.append(a)
    mass = (entry.get("sequence") or {}).get("molWeight")
    n_pdb = sum(1 for x in entry.get("uniProtKBCrossReferences", []) or [] if x.get("database") == "PDB")
    return {
        "aliases": "; ".join(aliases[:12]), "mass_kda": round(mass / 1000, 1) if mass else None, "n_pdb": n_pdb,
        "accession": entry.get("primaryAccession"),
        "gene": gene,
        "name": name,
        "length": (entry.get("sequence") or {}).get("length"),
        "segments": segments,
        "n_tm": n_tm,
        "gpi": gpi, "gpi_site": gpi_site,
        "signal": signal,
        "function": texts["function"], "location": texts["location"], "tissue": texts["tissue"],
        "ptm": texts["ptm"], "n_glyc": n_glyc, "n_ss": n_ss, "n_isoforms": isoforms, "shedding": shed,
    }


def _comment_texts(entry):
    """Collect the free-text UniProt comments we show to the user."""
    out = {"function": [], "location": [], "tissue": [], "ptm": []}
    for c in entry.get("comments", []):
        kind = c.get("commentType")
        if kind == "FUNCTION":
            out["function"] += [t.get("value", "") for t in c.get("texts", [])]
        elif kind == "SUBCELLULAR LOCATION":
            out["location"] += [loc.get("location", {}).get("value", "") for loc in c.get("subcellularLocations", [])]
        elif kind == "TISSUE SPECIFICITY":
            out["tissue"] += [t.get("value", "") for t in c.get("texts", [])]
        elif kind == "PTM":
            out["ptm"] += [t.get("value", "") for t in c.get("texts", [])]
    return {k: " ".join(x for x in v if x).strip() for k, v in out.items()}


def longest_extracellular(segments):
    lengths = [e - s + 1 for k, s, e, _ in segments if k == "extracellular" and s and e]
    return max(lengths) if lengths else 0


def classify(parsed):
    """Estimate surface accessibility from topology. Returns (target_class, accessibility).

    Rules (heuristic, shown to the user as an estimate):
      GPI-anchored                         -> fully outside the cell: yes
      one transmembrane helix              -> yes if the outside segment is 30 aa or longer, else partial
      several transmembrane helices        -> yes if an outside segment is 60 aa or longer, else partial
      membrane protein, no outside segment -> unknown (topology not annotated)
      no membrane anchor                   -> not a surface protein
    """
    longest = longest_extracellular(parsed["segments"])
    if parsed["gpi"]:
        return "gpi_anchored", "yes"
    if parsed["n_tm"] == 0:
        raise NotSurfaceProtein(
            f"UniProt lists no transmembrane helix or GPI anchor for {parsed['gene']}, "
            "so it cannot be used as a live-cell surface marker.")
    if longest == 0:
        return ("single_pass_tm" if parsed["n_tm"] == 1 else "multi_pass_tm"), "unknown"
    if parsed["n_tm"] == 1:
        return "single_pass_tm", ("yes" if longest >= 30 else "partial")
    return "multi_pass_tm", ("yes" if longest >= 60 else "partial")


def add_or_refresh_marker(con, symbol, species="human", role=None, note=None, basis="hypothesis"):
    """Fetch topology and store it. Creates the target if needed. Returns (ok, message)."""
    symbol = symbol.strip()
    if not symbol:
        return False, "Enter a gene symbol."
    try:
        entry = fetch_uniprot(symbol, species)
    except SourceUnavailable:
        return False, "Could not reach UniProt. Check your internet connection and try again."
    if entry is None:
        return False, f"No reviewed UniProt entry found for {symbol} ({species}). Check the spelling."
    parsed = parse_entry(entry)
    try:
        tclass, access = classify(parsed)
    except NotSurfaceProtein as exc:
        return False, str(exc)
    verified = int((parsed["gene"] or "").upper() == symbol.upper())
    row = con.execute("SELECT target_id, live_cell_accessible FROM target WHERE symbol = ? AND species = ?",
                      (symbol, species)).fetchone()
    if row is None:
        con.execute(
            "INSERT INTO target (symbol, protein_name, species, uniprot_acc, uniprot_verified, target_class, "
            "live_cell_accessible, accessibility_note) VALUES (?,?,?,?,?,?,?,?)",
            (symbol, parsed["name"], species, parsed["accession"], verified, tclass, access,
             "Estimated automatically from UniProt topology."))
        target_id = con.execute("SELECT last_insert_rowid()").fetchone()[0]
        created = True
    else:
        target_id, created = row["target_id"], False
        # Keep hand-entered accessibility judgments; only fill in what is still unknown.
        con.execute("UPDATE target SET uniprot_acc = ?, uniprot_verified = ?, protein_name = COALESCE(protein_name, ?), "
                    "target_class = CASE WHEN live_cell_accessible = 'unknown' THEN ? ELSE target_class END, "
                    "live_cell_accessible = CASE WHEN live_cell_accessible = 'unknown' THEN ? ELSE live_cell_accessible END "
                    "WHERE target_id = ?",
                    (parsed["accession"], verified, parsed["name"], tclass, access, target_id))
    con.execute("DELETE FROM target_region WHERE target_id = ? AND source LIKE 'UniProt%'", (target_id,))
    today = date.today().isoformat()
    segs = parsed["segments"]
    if parsed["gpi"] and parsed["length"]:
        site = parsed["gpi_site"] or parsed["length"]
        segs = segs + [("extracellular", 1, site, "GPI-anchored: the mature protein sits outside the cell")]
        if site < parsed["length"]:
            segs.append(("other", site + 1, parsed["length"], "C-terminal GPI signal, removed in the mature protein"))
    for kind, start, end, seg_note in segs:
        con.execute("INSERT INTO target_region (target_id, region_type, start_aa, end_aa, source, retrieved_on, note) "
                    "VALUES (?,?,?,?,?,?,?)", (target_id, kind, start, end, f"UniProt {parsed['accession']}", today, seg_note))
    store_info(con, target_id, parsed)
    if role:
        con.execute("INSERT OR IGNORE INTO cpec_context (target_id, context, evidence_basis, note) VALUES (?,?,?,?)",
                    (target_id, role, basis, note or "Added by user. No source recorded."))
    con.commit()
    longest = longest_extracellular(segs)
    verb = "Added" if created else "Updated"
    return True, f"{verb} {symbol}: {tclass.replace('_', ' ')}, longest outside segment {longest} aa, estimated access {access}."


# ---------- Europe PMC ----------
def search_literature(con, target_id, symbol, protein_name=None, limit=10):
    """Find papers mentioning the marker with choroid plexus and sorting or surface terms."""
    names = f'"{symbol}"' + (f' OR "{protein_name}"' if protein_name else "")
    query = (f'("choroid plexus") AND ({names}) AND '
             '("flow cytometry" OR FACS OR "cell sorting" OR "surface marker" OR "cell surface")')
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + urllib.parse.quote(query)
           + f"&format=json&pageSize={limit}")
    results = _get_json(url).get("resultList", {}).get("result", [])
    today = date.today().isoformat()
    n = 0
    for r in results:
        sid = r.get("id") or r.get("pmid")
        if not sid:
            continue
        con.execute("INSERT OR IGNORE INTO literature_hit (target_id, source_id, title, year, journal, doi, fetched_on) "
                    "VALUES (?,?,?,?,?,?,?)",
                    (target_id, str(sid), r.get("title"), r.get("pubYear"), r.get("journalTitle"), r.get("doi"), today))
        n += 1
    con.commit()
    return n


def refresh_all(con, progress=None):
    """Refresh topology and papers for every membrane target. Returns a list of message strings."""
    rows = con.execute("SELECT target_id, symbol, species, protein_name FROM target "
                       "WHERE target_class IN ('single_pass_tm','multi_pass_tm','gpi_anchored')").fetchall()
    messages = []
    for i, r in enumerate(rows):
        if progress:
            progress(i / max(len(rows), 1), f"Updating {r['symbol']}")
        ok, msg = add_or_refresh_marker(con, r["symbol"], r["species"])
        if not ok and "Could not reach" in msg:
            return [msg]
        if not ok:
            messages.append(f"{r['symbol']}: {msg}")
            continue
        try:
            n = search_literature(con, r["target_id"], r["symbol"], r["protein_name"])
            try:
                got_hpa = refresh_hpa(con, r["target_id"], r["symbol"]) if r["species"] == "human" else False
            except SourceUnavailable:
                got_hpa = False
            messages.append(f"{r['symbol']}: facts updated, {n} papers found"
                            + (", Protein Atlas data added." if got_hpa else "."))
        except SourceUnavailable:
            messages.append(f"{r['symbol']}: topology updated, could not reach Europe PMC for papers.")
        time.sleep(0.2)
    if progress:
        progress(1.0, "Done")
    return messages


# ---------- Marker facts ----------
def store_info(con, target_id, parsed, hpa=None):
    """Save the descriptive facts for a marker. Keeps earlier Protein Atlas values if none are supplied."""
    old = con.execute("SELECT hpa_ensembl, hpa_cp_ntpm, hpa_location FROM target_info WHERE target_id = ?",
                      (target_id,)).fetchone()
    hpa = hpa or (dict(old) if old else {})
    con.execute(
        "INSERT INTO target_info (target_id, function_text, location_text, tissue_text, ptm_text, n_glycosylation, "
        "n_disulfide, n_isoforms, shedding_hint, hpa_ensembl, hpa_cp_ntpm, hpa_location, fetched_on, "
        "aliases, mass_kda, n_pdb, length_aa, n_tm, has_signal) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(target_id) DO UPDATE SET "
        "aliases=excluded.aliases, mass_kda=excluded.mass_kda, n_pdb=excluded.n_pdb, length_aa=excluded.length_aa, "
        "n_tm=excluded.n_tm, has_signal=excluded.has_signal, "
        "function_text=excluded.function_text, location_text=excluded.location_text, tissue_text=excluded.tissue_text, "
        "ptm_text=excluded.ptm_text, n_glycosylation=excluded.n_glycosylation, n_disulfide=excluded.n_disulfide, "
        "n_isoforms=excluded.n_isoforms, shedding_hint=excluded.shedding_hint, hpa_ensembl=excluded.hpa_ensembl, "
        "hpa_cp_ntpm=excluded.hpa_cp_ntpm, hpa_location=excluded.hpa_location, fetched_on=excluded.fetched_on",
        (target_id, parsed["function"][:900], parsed["location"][:300], parsed["tissue"][:900], parsed["ptm"][:900],
         parsed["n_glyc"], parsed["n_ss"], parsed["n_isoforms"], int(parsed["shedding"]),
         hpa.get("hpa_ensembl"), hpa.get("hpa_cp_ntpm"), hpa.get("hpa_location"), date.today().isoformat(),
         parsed.get("aliases"), parsed.get("mass_kda"), parsed.get("n_pdb"), parsed.get("length"), parsed.get("n_tm"),
         int(bool(parsed.get("signal")))))
    con.commit()


def _pick(record, *needles):
    """Value of the first key containing every needle (case-insensitive). Column titles vary, so match loosely."""
    for key, value in record.items():
        low = key.lower()
        if all(n in low for n in needles):
            return value
    return None


def fetch_hpa(symbol):
    """Human Protein Atlas: Ensembl ID, choroid plexus RNA level, main subcellular location. Returns a dict or None."""
    url = ("https://www.proteinatlas.org/api/search_download.php?search=" + urllib.parse.quote(symbol)
           + "&format=json&columns=g,eg,t_RNA_choroid_plexus,scml,scal&compress=no")
    rows = _get_json(url)
    if not isinstance(rows, list):
        return None
    for rec in rows:
        gene = _pick(rec, "gene") if "Gene" not in rec else rec["Gene"]
        if str(gene or "").upper() == symbol.upper():
            cp = _pick(rec, "choroid")
            loc = _pick(rec, "main location") or _pick(rec, "subcellular")
            ens = rec.get("Ensembl") or _pick(rec, "ensembl")
            return {"hpa_ensembl": ens, "hpa_cp_ntpm": None if cp is None else str(cp),
                    "hpa_location": ", ".join(loc) if isinstance(loc, list) else (str(loc) if loc else None)}
    return None


def refresh_hpa(con, target_id, symbol):
    """Add Protein Atlas facts to a marker that already has UniProt facts. Returns True if anything was found."""
    hpa = fetch_hpa(symbol)
    if not hpa:
        return False
    con.execute("UPDATE target_info SET hpa_ensembl = ?, hpa_cp_ntpm = ?, hpa_location = ? WHERE target_id = ?",
                (hpa["hpa_ensembl"], hpa["hpa_cp_ntpm"], hpa["hpa_location"], target_id))
    con.commit()
    return True


# ---------- Antibody Registry (vendor lookup) ----------
FLOW_WORDS = re.compile(r"\b(flow|facs|fc)\b", re.I)


def registry_search(symbol, species="human", size=60):
    """Search the Antibody Registry (open, no key). Returns normalized records for this target and species.
    The registry lists vendor, catalog number, RRID, clone, conjugate, listed applications and citations."""
    url = ("https://www.antibodyregistry.org/api/fts-antibodies?q=" + urllib.parse.quote(symbol)
           + f"&size={size}&page=1")
    items = _get_json(url).get("items", [])
    out = []
    for it in items:
        target = str(it.get("abTarget") or "")
        name = str(it.get("abName") or "")
        if target.upper() != symbol.upper() and symbol.upper() not in name.upper().replace("-", "").replace(" ", ""):
            continue
        sp = str(it.get("targetSpecies") or "")
        if sp and species.lower() not in sp.lower():
            continue
        apps = it.get("applications")
        apps = ", ".join(apps) if isinstance(apps, list) else str(apps or "")
        ab_id = it.get("abId") or it.get("ab_id")
        rrid = f"AB_{ab_id}" if ab_id else (it.get("rrid") or "")
        out.append({
            "rrid": rrid, "name": name, "vendor": it.get("vendorName") or "", "catalog": str(it.get("catalogNum") or ""),
            "clone": it.get("cloneId") or "", "clonality": (it.get("clonality") or "").lower(),
            "conjugate": it.get("productConjugate") or "", "isotype": it.get("productIsotype") or "",
            "host": it.get("sourceOrganism") or "", "epitope": it.get("epitope") or "", "species": sp,
            "applications": apps, "flow_listed": bool(FLOW_WORDS.search(apps)),
            "citations": it.get("numOfCitation") or 0, "vendor_url": it.get("vendorUrl") or "",
            "url": it.get("url") or "",
        })
    return out


# ---------- Discover more markers ----------
def discover_cp_markers(known_symbols=(), min_ntpm=10.0, limit=60):
    """Plasma membrane proteins (Human Protein Atlas) ranked by choroid plexus RNA level.
    RNA level is a lead, not proof of surface protein or selectivity: many of these are also made in other tissues."""
    query = urllib.parse.quote("subcell_location:Plasma membrane")
    url = ("https://www.proteinatlas.org/api/search_download.php?search=" + query
           + "&format=json&columns=g,eg,up,t_RNA_choroid_plexus,scml&compress=no")
    rows = _get_json(url, timeout=60)
    if not isinstance(rows, list):
        return []
    known = {k.upper() for k in known_symbols}
    out = []
    for rec in rows:
        gene = rec.get("Gene") or _pick(rec, "gene")
        raw = _pick(rec, "choroid")
        try:
            level = float(raw)
        except (TypeError, ValueError):
            continue
        if not gene or str(gene).upper() in known or level < min_ntpm:
            continue
        loc = _pick(rec, "main location")
        out.append({"symbol": str(gene), "ensembl": rec.get("Ensembl") or _pick(rec, "ensembl"),
                    "cp_ntpm": level, "location": ", ".join(loc) if isinstance(loc, list) else (loc or "")})
    out.sort(key=lambda r: -r["cp_ntpm"])
    return out[:limit]
