"""Tests for the automatic parts, using mock UniProt and Europe PMC responses (no internet needed).
Run:  python test_logic.py
The mock topology numbers are invented for testing and are not real protein annotations."""
import shutil
import sys
import tempfile
from pathlib import Path

import autofill
import db
import ranking

tmp = Path(tempfile.mkdtemp())
for f in Path(__file__).parent.iterdir():
    if f.is_file() and f.suffix in (".sql", ".md", ".py", ".txt", ".cff", ".json", ".sh", ".bat") or f.name in ("LICENSE", "VERSION"):
        shutil.copy(f, tmp / f.name)
db.HERE = tmp
db.DB = tmp / "test.sqlite"


def feat(kind, start, end, desc=""):
    return {"type": kind, "description": desc, "location": {"start": {"value": start}, "end": {"value": end}}}


def entry(gene, acc, feats, length=300, name="Test protein"):
    return {"primaryAccession": acc, "genes": [{"geneName": {"value": gene}}], "sequence": {"length": length},
            "proteinDescription": {"recommendedName": {"fullName": {"value": name}}}, "features": feats}


MOCK = {
    "FOLR1": entry("FOLR1", "P15328", [feat("Lipidation", 233, 233, "GPI-anchor amidated serine"),
                                       feat("Glycosylation", 69, 69, "N-linked (GlcNAc...) asparagine"),
                                       feat("Glycosylation", 161, 161, "N-linked (GlcNAc...) asparagine"),
                                       feat("Glycosylation", 201, 201, "N-linked (GlcNAc...) asparagine"),
                                       feat("Disulfide bond", 28, 55, ""), feat("Disulfide bond", 48, 96, "")], 257,
                   "Folate receptor alpha"),
    "NEWSURF": entry("NEWSURF", "Q00001", [feat("Topological domain", 1, 200, "Extracellular"),
                                           feat("Transmembrane", 201, 221, "Helical"),
                                           feat("Topological domain", 222, 260, "Cytoplasmic")], 260),
    "SHORTLOOP": entry("SHORTLOOP", "Q00002", [feat("Topological domain", 1, 10, "Cytoplasmic"),
                                               feat("Transmembrane", 11, 30, "Helical"),
                                               feat("Topological domain", 31, 50, "Extracellular"),
                                               feat("Transmembrane", 51, 70, "Helical"),
                                               feat("Topological domain", 71, 120, "Cytoplasmic")], 120),
    "DISC1": entry("DISC1", "Q00004", [feat("Topological domain", 1, 120, "Extracellular"), feat("Transmembrane", 121, 141, "Helical"),
                                       feat("Topological domain", 142, 180, "Cytoplasmic")], 180),
    "SECRETED": entry("SECRETED", "Q00003", [feat("Signal", 1, 20)], 150),
}
autofill.fetch_uniprot = lambda sym, sp="human": MOCK.get(sym)

MOCK["FOLR1"]["comments"] = [
    {"commentType": "FUNCTION", "texts": [{"value": "Binds folate and mediates its uptake."}]},
    {"commentType": "SUBCELLULAR LOCATION", "subcellularLocations": [{"location": {"value": "Cell membrane"}}]},
    {"commentType": "TISSUE SPECIFICITY", "texts": [{"value": "Choroid plexus and kidney."}]},
    {"commentType": "PTM", "texts": [{"value": "Soluble form is released by proteolytic cleavage and shedding."}]},
    {"commentType": "ALTERNATIVE PRODUCTS", "isoforms": [{"name": "1"}, {"name": "2"}]}]
REGISTRY = {"items": [
    {"abId": 1001, "abName": "Folate Receptor alpha Antibody (FITC)", "abTarget": "FOLR1", "clonality": "monoclonal",
     "cloneId": "LK26", "productConjugate": "FITC", "productIsotype": "IgG1", "sourceOrganism": "mouse",
     "epitope": "Extracellular domain, residues 25-233", "targetSpecies": "Human", "applications": ["Flow Cytometry", "IHC"],
     "vendorName": "VendorA", "catalogNum": "A-100", "numOfCitation": 6, "vendorUrl": "https://example.com/a"},
    {"abId": 1002, "abName": "FOLR1 antibody C-terminal cytoplasmic", "abTarget": "FOLR1", "clonality": "polyclonal",
     "epitope": "Cytoplasmic tail residues 235-257", "targetSpecies": "Human", "applications": "WB",
     "vendorName": "VendorB", "catalogNum": "B-200", "numOfCitation": 0},
    {"abId": 1003, "abName": "Unrelated antibody", "abTarget": "GAPDH", "targetSpecies": "Human"},
    {"abId": 1004, "abName": "FOLR1 mouse only", "abTarget": "FOLR1", "targetSpecies": "Mouse", "vendorName": "VendorC"}]}
HPA = [{"Gene": "FOLR1", "Ensembl": "ENSG00000110195", "Tissue RNA - choroid plexus [nTPM]": 812.4,
        "Subcellular main location": ["Plasma membrane"]}, {"Gene": "FOLR2", "Ensembl": "ENSG00000165457"},
       {"Gene": "NEWGENE1", "Ensembl": "ENSG1", "Tissue RNA - choroid plexus [nTPM]": 500.0, "Subcellular main location": ["Plasma membrane"]},
       {"Gene": "MIDGENE", "Ensembl": "ENSG2", "Tissue RNA - choroid plexus [nTPM]": 120.5, "Subcellular main location": ["Plasma membrane"]},
       {"Gene": "LOWGENE", "Ensembl": "ENSG3", "Tissue RNA - choroid plexus [nTPM]": 2.0},
       {"Gene": "NOVALUE", "Ensembl": "ENSG4", "Tissue RNA - choroid plexus [nTPM]": None}]
def fake_json(url, timeout=25):
    if "antibodyregistry" in url:
        return REGISTRY
    if "proteinatlas" in url:
        return HPA
    return {"resultList": {"result": [{"id": "12345", "pubYear": "2024", "journalTitle": "J Test",
                                       "title": "Sorting choroid plexus cells", "doi": "10.1/x"}]}}
autofill._get_json = fake_json

failures = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        failures.append(name)


db.ensure_db()
con = db.connect()

# classification rules
check("GPI -> yes", autofill.classify(autofill.parse_entry(MOCK["FOLR1"])) == ("gpi_anchored", "yes"))
check("single pass long outside -> yes", autofill.classify(autofill.parse_entry(MOCK["NEWSURF"])) == ("single_pass_tm", "yes"))
check("multi pass short loops -> partial", autofill.classify(autofill.parse_entry(MOCK["SHORTLOOP"])) == ("multi_pass_tm", "partial"))
try:
    autofill.classify(autofill.parse_entry(MOCK["SECRETED"]))
    check("secreted protein rejected", False)
except autofill.NotSurfaceProtein:
    check("secreted protein rejected", True)

# adding a new marker, and refreshing a seeded one without losing curated judgments
ok, msg = autofill.add_or_refresh_marker(con, "NEWSURF", "human", "CPEC_identity")
check("new marker added", ok and con.execute("SELECT COUNT(*) FROM target WHERE symbol='NEWSURF'").fetchone()[0] == 1)
check("new marker has role", con.execute("SELECT COUNT(*) FROM cpec_context c JOIN target t USING(target_id) WHERE t.symbol='NEWSURF'").fetchone()[0] == 1)
ok, msg = autofill.add_or_refresh_marker(con, "SECRETED", "human")
check("secreted marker refused with reason", (not ok) and "cannot be used" in msg)
ok, msg = autofill.add_or_refresh_marker(con, "NOSUCHGENE", "human")
check("unknown gene reported", (not ok) and "No reviewed UniProt entry" in msg)
before = con.execute("SELECT live_cell_accessible, accessibility_note FROM target WHERE symbol='FOLR1'").fetchone()
ok, _ = autofill.add_or_refresh_marker(con, "FOLR1", "human")
after = con.execute("SELECT live_cell_accessible, accessibility_note, uniprot_verified FROM target WHERE symbol='FOLR1'").fetchone()
check("refresh keeps curated note and marks verified", ok and before["accessibility_note"] == after["accessibility_note"] and after["uniprot_verified"] == 1)
n = autofill.search_literature(con, con.execute("SELECT target_id FROM target WHERE symbol='FOLR1'").fetchone()[0], "FOLR1")
check("literature stored", n == 1 and con.execute("SELECT COUNT(*) FROM literature_hit").fetchone()[0] == 1)
autofill.search_literature(con, con.execute("SELECT target_id FROM target WHERE symbol='FOLR1'").fetchone()[0], "FOLR1")
check("literature not duplicated", con.execute("SELECT COUNT(*) FROM literature_hit").fetchone()[0] == 1)

# immunogen check
tid = con.execute("SELECT target_id FROM target WHERE symbol='NEWSURF'").fetchone()[0]
check("immunogen outside -> extracellular", ranking.check_immunogen(con, tid, 20, 80)[0] == "extracellular")
check("immunogen inside -> intracellular", ranking.check_immunogen(con, tid, 230, 255)[0] == "intracellular")
check("immunogen across -> mixed", ranking.check_immunogen(con, tid, 150, 230)[0] == "mixed")
check("no residues -> unknown", ranking.check_immunogen(con, tid, 0, 0)[0] == "unknown")

# ranking
ranked = ranking.rank_targets(con)
syms = [r["symbol"] for r in ranked]
check("ranking returns items", len(ranked) > 5)
check("ubiquitous pump ranks below FOLR1", syms.index("ATP1A1") > syms.index("FOLR1"))
check("nothing proven in a fresh database", not any(r["proven"] for r in ranked))
check("strategy has keep and remove", {r["Action"] for r in ranking.suggest_sort_strategy(ranked)} == {"Keep", "Remove"})
rf = [r for r in ranked if r["symbol"] == "FOLR1"][0]
check("papers attached to FOLR1", len(rf["papers"]) == 1)

# a logged live test makes a marker 'proven'
aid = con.execute("INSERT INTO antibody (target_id, clone, curation_status) VALUES (?,?,?)", (tid, "X1", "needs_curation")).lastrowid
con.execute("INSERT INTO validation (antibody_id, application, nonpermeabilized, context_matches_cpec, outcome) VALUES (?,?,?,?,?)",
            (aid, "live_flow", 1, 1, "works"))
con.commit()
check("proven after a passed live test", [r for r in ranking.rank_targets(con) if r["symbol"] == "NEWSURF"][0]["proven"])
# inside-the-cell antibodies are hidden
con.execute("INSERT INTO epitope (target_id, epitope_type, region_type, start_aa, end_aa) VALUES (?,?,?,?,?)", (tid, "linear", "cytoplasmic", 230, 255))
eid = con.execute("SELECT last_insert_rowid()").fetchone()[0]
con.execute("INSERT INTO antibody (target_id, epitope_id, clone, curation_status) VALUES (?,?,?,?)", (tid, eid, "INSIDE", "needs_curation"))
con.commit()
labels = [a["label"] for a in [r for r in ranking.rank_targets(con) if r["symbol"] == "NEWSURF"][0]["antibodies"]]
check("inside-the-cell antibody hidden", "INSIDE" not in labels and "X1" in labels)

# ---------- marker facts ----------
folr = con.execute("SELECT target_id FROM target WHERE symbol='FOLR1'").fetchone()[0]
info = con.execute("SELECT * FROM target_info WHERE target_id = ?", (folr,)).fetchone()
check("facts stored", info is not None and info["n_glycosylation"] == 3 and info["n_disulfide"] == 2 and info["n_isoforms"] == 2)
check("function and tissue text stored", "folate" in info["function_text"].lower() and "choroid" in info["tissue_text"].lower())
check("shedding hint detected", info["shedding_hint"] == 1)
check("protein atlas facts added", autofill.refresh_hpa(con, folr, "FOLR1"))
info = con.execute("SELECT * FROM target_info WHERE target_id = ?", (folr,)).fetchone()
check("protein atlas values stored", info["hpa_ensembl"] == "ENSG00000110195" and info["hpa_cp_ntpm"] == "812.4"
      and "Plasma" in info["hpa_location"])
ok, _ = autofill.add_or_refresh_marker(con, "FOLR1", "human")
info = con.execute("SELECT * FROM target_info WHERE target_id = ?", (folr,)).fetchone()
check("refresh keeps protein atlas values", info["hpa_cp_ntpm"] == "812.4")
t = con.execute("SELECT * FROM target WHERE symbol='FOLR1'").fetchone()
impl = ranking.marker_implications(t, info, 257, "apical")
check("implications mention glycans, shedding, polarity, atlas",
      all(any(k in x for x in impl) for k in ("glycosylation", "shedding", "apical", "Human Protein Atlas")))

# ---------- vendor lookup ----------
hits = autofill.registry_search("FOLR1", "human")
check("registry keeps only matching human FOLR1 records", [h["rrid"] for h in hits] == ["AB_1001", "AB_1002"])
check("registry fields parsed", hits[0]["vendor"] == "VendorA" and hits[0]["catalog"] == "A-100" and hits[0]["clone"] == "LK26"
      and hits[0]["flow_listed"] and not hits[1]["flow_listed"] and hits[0]["citations"] == 6)
v0 = ranking.immunogen_from_text(con, folr, hits[0]["epitope"])
v1 = ranking.immunogen_from_text(con, folr, hits[1]["epitope"])
check("immunogen from residues: outside", v0[0] == "extracellular")
check("immunogen in GPI propeptide detected as removed", v1[0] == "removed")
check("keyword fallback", ranking.immunogen_from_text(con, folr, "ectodomain fragment")[0] == "extracellular")
check("terminus alone stays unknown", ranking.immunogen_from_text(con, folr, "C-terminal region")[0] == "unknown")
f0 = ranking.antibody_flags({**hits[0], "azide_free": "unknown"}, v0[0], v0[1])
f1 = ranking.antibody_flags({**hits[1], "azide_free": "unknown"}, v1[0], v1[1])
check("good antibody has no problems", not any(x[0] == "problem" for x in f0))
check("removed-part antibody flagged as problem", any(x[0] == "problem" and "removed" in x[1].lower() for x in f1))
f_in = ranking.antibody_flags({"conjugate": "FITC"}, "intracellular", "x")
check("inside-the-cell antibody flagged as problem", any(x[0] == "problem" and "inside" in x[1].lower() for x in f_in))
check("azide 'no' flagged", any(x[0] == "problem" and "azide" in x[1].lower()
      for x in ranking.antibody_flags({**hits[0], "azide_free": "no"}, "extracellular", "")))
check("scores order good above bad", ranking.antibody_points({**hits[0], "azide_free": "unknown"}, v0[0])
      > ranking.antibody_points({**hits[1], "azide_free": "unknown"}, v1[0]))

# ---------- staining plan ----------
plan = ranking.test_plan("FOLR1", "Keep: marks CPECs", {"label": "LK26 VendorA A-100", "isotype": "IgG1", "conjugate": "FITC"},
                         "iPSC-derived CPECs day 30", True, "apical", False)
check("plan has controls and key steps", all(k in plan for k in ("Unstained", "Isotype", "Viability", "Accutase", "azide", "Do NOT fix")))
check("plan has marker cautions", "shed" in plan and "apical" in plan)
check("contaminant plan uses CPEC-directed negative", "CPEC-directed cells" in
      ranking.test_plan("PODXL", "Remove", None, "x", False, None, True))

# ---------- order list storage ----------
# ---------- discovery ----------
found = autofill.discover_cp_markers(known_symbols=["FOLR1", "FOLR2"])
check("discovery ranks by choroid plexus RNA", [f["symbol"] for f in found] == ["NEWGENE1", "MIDGENE"])
check("discovery skips known, low and empty", all(f["symbol"] not in ("FOLR1", "LOWGENE", "NOVALUE") for f in found))
ok, _ = autofill.add_or_refresh_marker(con, "DISC1", "human", "CPEC_identity", note="RNA lead only.", basis="hypothesis")
row = con.execute("SELECT c.note, c.evidence_basis FROM cpec_context c JOIN target t USING(target_id) WHERE t.symbol='DISC1'").fetchone()
check("discovered marker keeps its note and basis", ok and row["note"] == "RNA lead only." and row["evidence_basis"] == "hypothesis")

# ---------- evidence labels ----------
ranked = {r["symbol"]: r for r in ranking.rank_targets(con)}
check("species carried on each marker", all(r["species"] in ("human", "mouse") for r in ranked.values()))
check("source-supported marker labeled as such", ranked["AQP1"]["evidence"]["level"] == "search_verified"
      and "not a live-cell test" in ranked["AQP1"]["evidence"]["text"])
check("knowledge-only marker labeled as such", ranked["ATP1A1"]["evidence"]["level"] == "prior_knowledge")
check("hypothesis marker labeled as such", ranked["DISC1"]["evidence"]["level"] == "hypothesis")
check("reported live sort elsewhere is named", ranked["CD99"]["evidence"]["level"] == "reported" and "CD99_DA_2022" in ranked["CD99"]["evidence"]["text"])
check("reported sort does not count as proven in CPECs", not ranked["CD99"]["proven"])
check("own test summarized", ranked["NEWSURF"]["evidence"]["level"] == "own" and "1 worked" in ranked["NEWSURF"]["evidence"]["text"])
rows = ranking.suggest_sort_strategy(list(ranked.values()))
check("strategy shows evidence text, never a default 'not yet'", all("Evidence" in r and "not yet" not in str(r).lower() for r in rows))
# a failed own test is reported as failed, not as untested
aid2 = con.execute("INSERT INTO antibody (target_id, clone, curation_status) VALUES (?,?,?)", (folr, "F1", "needs_curation")).lastrowid
con.execute("INSERT INTO validation (antibody_id, application, nonpermeabilized, context_matches_cpec, outcome) VALUES (?,?,?,?,?)", (aid2, "live_flow", 1, 1, "fails"))
con.commit()
ev = [r for r in ranking.rank_targets(con) if r["symbol"] == "FOLR1"][0]["evidence"]
check("own failure reported", ev["level"] == "own" and "1 did not work" in ev["text"])


# ---------- readable names and informative fields ----------
import labels
check("Title Case rule", [labels.tc(x) for x in ("update from public databases", "add to order list", "what is it for?",
      "live-cell staining plan: FOLR1", "gene symbol (for example FOLR1)", "choroid plexus RNA (nTPM)")] ==
      ["Update From Public Databases", "Add to Order List", "What Is It For?", "Live-Cell Staining Plan: FOLR1",
       "Gene Symbol (For Example FOLR1)", "Choroid Plexus RNA (nTPM)"])
check("status codes read as labels", labels.pretty(labels.STATUS, "seed_unverified") == "Clone Name Only"
      and labels.pretty(labels.STATUS, "lab_verified") == "Passed in Your Lab")
check("unknown codes fall back to Title Case", labels.pretty(labels.REGION, "some_new_code") == "Some New Code")
rk = ranking.rank_targets(con)
check("score breakdown sums to the score", all(abs(sum(b["Points"] for b in r["breakdown"]) - r["score"]) < 0.11 for r in rk))
check("ranked items carry readable fields", all(r["access"] and r["specificity"] and r["klass"] and r["best_status"] for r in rk))
check("role labels are Title Case", all(labels.tc(r["role_label"]) == r["role_label"] for r in rk))
check("flag labels are Title Case", all(labels.tc(f[1]) == f[1] for f in
      ranking.antibody_flags({"conjugate": "FITC", "clonality": "monoclonal", "citations": 3, "flow_listed": True}, "extracellular", "x")
      + ranking.antibody_flags({}, "removed", "x") + ranking.antibody_flags({"azide_free": "no"}, "unknown", "x")))
print()
print(f"{len(failures)} failure(s)")
sys.exit(1 if failures else 0)

# panel designer
import panels
_rows, _w = panels.check_panel([{"symbol": "A", "role": "keep", "conjugate": "PE"}, {"symbol": "B", "role": "remove", "conjugate": "PE"},
                                {"symbol": "C", "role": "remove", "conjugate": "", "azide_free": "no", "verdict": "intracellular"}])
_t = " ".join(_w)
print(("PASS" if all(x in _t for x in ("Channel clash", "unconjugated", "azide", "not outside")) else "FAIL") + " panel checks")

# extra marker facts
import autofill as _af
_p = _af.parse_entry({"primaryAccession": "Q0", "genes": [{"geneName": {"value": "EPCAM"}, "synonyms": [{"value": "TACSTD1"}]}],
                      "sequence": {"length": 314, "molWeight": 34932},
                      "proteinDescription": {"recommendedName": {"fullName": {"value": "Epithelial cell adhesion molecule"}},
                                             "alternativeNames": [{"fullName": {"value": "CD326"}}]},
                      "uniProtKBCrossReferences": [{"database": "PDB"}, {"database": "PDB"}, {"database": "Ensembl"}], "features": []})
print(("PASS" if _p["aliases"] == "TACSTD1; CD326" and _p["mass_kda"] == 34.9 and _p["n_pdb"] == 2 else "FAIL") + " marker facts parsed")
