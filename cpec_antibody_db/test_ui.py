"""Screen tests: every visible label is Title Case, and the main flows work. Uses mock web responses and a
temporary copy of the database, so it needs no internet and changes nothing.  Run:  python test_ui.py"""
import re
import shutil
import sqlite3
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import autofill
import db
import labels

tmp = Path(tempfile.mkdtemp())
shutil.copy(HERE / "cpec_antibodies.sqlite", tmp / "ui.sqlite")
db.DB = tmp / "ui.sqlite"
db.ensure_db()


def feat(k, s, e, d=""):
    return {"type": k, "description": d, "location": {"start": {"value": s}, "end": {"value": e}}}


UNIPROT = {
    "NEWGENE1": {"primaryAccession": "Q1", "genes": [{"geneName": {"value": "NEWGENE1"}}], "sequence": {"length": 300},
                 "proteinDescription": {"recommendedName": {"fullName": {"value": "New gene one"}}},
                 "features": [feat("Topological domain", 1, 250, "Extracellular"), feat("Transmembrane", 251, 271, "Helical"),
                              feat("Topological domain", 272, 300, "Cytoplasmic")]},
    "MIDGENE": {"primaryAccession": "Q2", "genes": [{"geneName": {"value": "MIDGENE"}}], "sequence": {"length": 100},
                "proteinDescription": {"recommendedName": {"fullName": {"value": "Secreted-like"}}}, "features": [feat("Signal", 1, 20)]}}
REGISTRY = {"items": [
    {"abId": 1001, "abName": "FOLR1 antibody FITC", "abTarget": "FOLR1", "clonality": "monoclonal", "cloneId": "LK26",
     "productConjugate": "FITC", "productIsotype": "IgG1", "epitope": "Extracellular domain, residues 25-233",
     "targetSpecies": "Human", "applications": ["Flow Cytometry"], "vendorName": "VendorA", "catalogNum": "A-100", "numOfCitation": 6},
    {"abId": 1002, "abName": "FOLR1 antibody", "abTarget": "FOLR1", "clonality": "polyclonal", "targetSpecies": "Human",
     "applications": "WB", "vendorName": "VendorB", "catalogNum": "B-200"}]}
HPA = [{"Gene": "NEWGENE1", "Ensembl": "E1", "Tissue RNA - choroid plexus [nTPM]": 500.0, "Subcellular main location": ["Plasma membrane"]},
       {"Gene": "MIDGENE", "Ensembl": "E2", "Tissue RNA - choroid plexus [nTPM]": 120.0, "Subcellular main location": ["Plasma membrane"]}]
autofill.fetch_uniprot = lambda s, sp="human": UNIPROT.get(s)
autofill._get_json = lambda url, timeout=25: (REGISTRY if "antibodyregistry" in url else HPA if "proteinatlas" in url
                                              else {"resultList": {"result": []}})

from streamlit.testing.v1 import AppTest

failures = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        failures.append(name)


def W(ws, label):
    return [w for w in ws if getattr(w, "label", None) == label][0]


def visible_labels(at):
    out = []
    for kind in ("button", "selectbox", "multiselect", "text_input", "text_area", "number_input", "checkbox", "radio",
                 "metric", "expander"):
        out += [(kind, w.label) for w in getattr(at, kind)]
    out += [("tab", t.label) for t in at.tabs]
    out += [("subheader", x.value) for x in at.subheader] + [("title", x.value) for x in at.title]
    out += [("download", e.proto.label) for e in at.get("download_button")]
    for d in at.dataframe:
        out += [("column", c) for c in d.value.columns]
    for m in at.markdown:
        for line in m.value.splitlines():
            if line.startswith("### "):
                out.append(("heading", line[4:]))
        out += [("bold", b.strip(" :.")) for b in re.findall(r"\*\*(.+?)\*\*", m.value)]
    return out


def lint(at, where):
    bad = [(k, l) for k, l in visible_labels(at) if l and labels.tc(l) != l]
    check(f"all labels are Title Case ({where})", not bad)
    for k, l in bad[:12]:
        print(f"      {k}: {l!r} should be {labels.tc(l)!r}")


at = AppTest.from_file(str(HERE / "app.py"), default_timeout=90).run()
check("loads without errors", not at.exception)
lint(at, "first screen")
c = sqlite3.connect(tmp / "ui.sqlite")

n_roles = c.execute("select count(distinct target_id) from cpec_context").fetchone()[0]
big = [d.value for d in at.dataframe if "Evidence on File" in d.value.columns][0]
check("every marker is listed, not five", len(big) == n_roles and len(big) > 5)
check("ranked table has the informative columns", all(x in big.columns for x in (
    "Estimated Live-Cell Access", "Longest Outside Segment (aa)", "Specificity", "Choroid Plexus RNA (nTPM)",
    "Papers Found", "Antibodies on File")))
check("values read as labels, not codes", not any("_" in str(v) for col in ("Role", "Estimated Live-Cell Access", "Specificity")
                                                 for v in big[col]))
check("summary metrics shown for the open marker", {"Estimated Live-Cell Access", "Longest Outside Segment", "Papers Found"}
      <= {m.label for m in at.metric})
check("score breakdown available", any(e.label == "Score Breakdown" for e in at.expander))
check("no cell-description question", not any("cells" in e.label.lower() and "my" in e.label.lower() for e in at.expander))
alltext = " ".join(m.value for m in at.markdown) + " ".join(str(d.value.to_dict()) for d in at.dataframe)
check("no default 'not yet' label", "not yet" not in alltext.lower())
strat = [d.value for d in at.dataframe if "Action" in d.value.columns][0]
check("strategy table has Title Case columns", list(strat.columns) == ["Marker", "Action", "Antibody", "Evidence"])
check("status codes shown as labels in records", all("_" not in str(v) for d in at.dataframe if "Status" in d.value.columns
                                                     for v in d.value["Status"]))

# log a result
W(at.selectbox, "Open a Marker to Test").select("FOLR1"); at.run()
tid = c.execute("select target_id from target where symbol='FOLR1'").fetchone()[0]
[r for r in at.radio if r.key == f"o_{tid}"][0].set_value("Worked")
[t for t in at.text_input if t.key == f"cl_{tid}"][0].input("H9-derived, day 30")
W(at.button, "Save Result").click(); at.run()
row = c.execute("select v.cell_context, v.species, v.outcome from validation v join antibody a using(antibody_id) "
                "where a.target_id=? and v.ref_id is null", (tid,)).fetchone()
check("result saved with the cell line and day", row and row[0] == "stem-cell-derived CPECs, H9-derived, day 30" and row[2] == "works")
check("stage updates when a result is saved", c.execute("select count(*) from antibody where stage='passed'").fetchone()[0] >= 1)
W(at.selectbox, "Open a Marker to Test").select("FOLR1"); at.run()
check("evidence label updates after a test", "Your tests: 1 worked" in " ".join(m.value for m in at.markdown))

# find antibodies
W(at.selectbox, "Marker").select("FOLR1")
W(at.button, "Search the Antibody Registry").click(); at.run()
reg = [d.value for d in at.dataframe if "Listed Applications" in d.value.columns]
check("registry results shown with Title Case values", bool(reg) and set(reg[0]["Clonality"]) <= {"Monoclonal", "Polyclonal", "Recombinant", "Unknown", ""})
lint(at, "after a registry search")
opts = W(at.multiselect, "Add Products to My Order List").options
W(at.multiselect, "Add Products to My Order List").select(opts[0])
W(at.button, "Add to Order List").click(); at.run()
check("product added to the order list", c.execute("select count(*) from antibody where rrid='AB_1001' and shortlisted=1").fetchone()[0] == 1)
order = [d.value for d in at.dataframe if "Immunogen Location" in d.value.columns and "Datasheet" in d.value.columns]
check("order list shows Title Case columns", bool(order))

# discovery
W(at.button, "Find Candidates").click(); at.run()
W(at.multiselect, "Add Markers").select("NEWGENE1").select("MIDGENE")
W(at.button, "Add Selected Markers").click(); at.run()
check("discovery adds surface proteins and refuses secreted ones",
      c.execute("select count(*) from target where symbol='NEWGENE1'").fetchone()[0] == 1
      and c.execute("select count(*) from target where symbol='MIDGENE'").fetchone()[0] == 0)
lint(at, "after discovery")

# new tabs and the records split
check("catalog loaded: many more markers", c.execute("select count(*) from target").fetchone()[0] >= 45)
check("no own tests listed before any are saved for literature-only markers",
      c.execute("select count(*) from validation where ref_id is null").fetchone()[0] == 1)
rec_tables = [d.value for d in at.dataframe if "Reference" in d.value.columns and "Result" in d.value.columns]
own_tables = [d.value for d in at.dataframe if "Source" in d.value.columns and "Result" in d.value.columns]
check("literature reports are in their own table", bool(rec_tables) and all(rec_tables[0]["Reference"].astype(str) != "Your Test"))
check("own tests table holds only own tests", bool(own_tables) and set(own_tables[0]["Source"]) == {"Your Test"})
check("pass counts split own and reported", any("Your Live-Cell Passes" in d.value.columns and "Reported in Literature" in d.value.columns for d in at.dataframe))
W(at.button, "Search the Registry for Every Marker").click(); at.run()
check("registry coverage stored", c.execute("select count(*) from registry_summary").fetchone()[0] > 0)
W(at.multiselect, "Keep Gate (Markers on the Cells You Want)").select(W(at.multiselect, "Keep Gate (Markers on the Cells You Want)").options[0]); at.run()
check("panel designer shows checks", any("Channel" in d.value.columns for d in at.dataframe))
W(at.text_input, "Run Date (YYYY-MM-DD)").input("2026-11-11")
W(at.button, "Save Run").click(); at.run()
check("experiment saved", c.execute("select count(*) from sort_experiment").fetchone()[0] == 1)
lint(at, "all new tabs")
W(at.selectbox, "Stage").select("Ordered")
W(at.text_input, "Lot Number").input("L123")
W(at.button, "Save Tracking").click(); at.run()
check("tracking saved", c.execute("select count(*) from antibody where lot='L123' and stage='ordered'").fetchone()[0] == 1)
W(at.text_input, "Dilution (For Example 1:100)").input("1:100")
at.number_input[[n.label for n in at.number_input].index("Positive Median")].set_value(1000.0)
at.number_input[[n.label for n in at.number_input].index("Negative Median")].set_value(100.0)
at.number_input[[n.label for n in at.number_input].index("Negative Standard Deviation")].set_value(50.0)
W(at.button, "Save Titration Point").click(); at.run()
check("titration saved", c.execute("select count(*) from titration").fetchone()[0] == 1)
check("staining index shown", any("Staining Index" in d.value.columns and 9.0 in d.value["Staining Index"].tolist() for d in at.dataframe))
lint(at, "tracker")
check("no errors at the end", not at.exception)
print(f"\n{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
