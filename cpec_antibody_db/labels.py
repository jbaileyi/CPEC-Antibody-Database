"""Display text: Title Case labels and plain-language values.

Stored values stay as they are in the database (for example seed_unverified). Only what the person sees changes.
Rebuilt from how app.py, ranking.py and the tests use it: tc() for Title Case, pretty(mapping, code) for values.
"""
import re

SMALL = {"a", "an", "and", "as", "at", "but", "by", "for", "if", "in", "of", "on", "or", "nor", "the", "to", "vs"}
ACRONYMS = {"CPEC": "CPEC", "CPECS": "CPECs", "RNA": "RNA", "RRID": "RRID", "DOI": "DOI", "UNIPROT": "UniProt",
            "PUBMED": "PubMed", "IPSC": "iPSC", "HESC": "hESC", "MACS": "MACS", "FACS": "FACS", "NTPM": "nTPM",
            "AA": "aa", "GPI": "GPI", "CSV": "CSV", "ID": "ID", "FITC": "FITC", "IGG": "IgG", "BIOLEGEND": "BioLegend",
            "CITEAB": "CiteAb", "BIOCOMPARE": "Biocompare", "FMO": "FMO", "ICC": "ICC", "IHC": "IHC", "OK": "OK"}
_TOKEN = re.compile(r"[A-Za-z0-9]+")


def _fix(tok, force):
    up = tok.upper()
    if up in ACRONYMS:
        return ACRONYMS[up]
    if any(ch.isdigit() for ch in tok):
        return tok
    if tok.isupper() and len(tok) > 1:
        return tok
    if any(c.isupper() for c in tok[1:]):
        return tok
    low = tok.lower()
    if not force and low in SMALL:
        return low
    return low.capitalize()


def tc(text):
    """Title Case for labels and headings. The first word, the last word, and any word after a colon or an opening
    parenthesis are always capitalized. Acronyms, gene symbols and mixed-case names are kept as written."""
    text = str(text)
    tokens = list(_TOKEN.finditer(text))
    out, last, force = [], 0, True
    for i, m in enumerate(tokens):
        between = text[last:m.start()]
        out.append(between)
        if re.search(r"[:(—]", between):
            force = True
        out.append(_fix(m.group(0), force or i == len(tokens) - 1))
        force = False
        last = m.end()
    out.append(text[last:])
    return "".join(out)


def pretty(mapping, code):
    """Readable form of a stored value. Unknown codes fall back to Title Case of the raw text."""
    if code is None or code == "":
        return ""
    return mapping.get(code, tc(str(code).replace("_", " ")))


STATUS = {"seed_unverified": "Clone Name Only", "needs_curation": "Needs Details",
          "vendor_verified": "Vendor Listed", "lab_verified": "Passed in Your Lab"}
ACCESS = {"yes": "Likely Accessible", "partial": "Partly Accessible", "unknown": "Unknown", "no": "Not Accessible"}
SPECIFICITY = {"selective": "Selective", "shared": "Shared With Other Cell Types",
               "ubiquitous": "Found Everywhere", "unknown": "Unknown"}
TARGET_CLASS = {"single_pass_tm": "Single-Pass Membrane Protein", "multi_pass_tm": "Multi-Pass Membrane Protein",
                "gpi_anchored": "GPI-Anchored Protein", "glycan_epitope": "Glycan Epitope",
                "glycoprotein_epitope": "Glycoprotein Epitope"}
CLONALITY = {"monoclonal": "Monoclonal", "polyclonal": "Polyclonal", "recombinant": "Recombinant", "unknown": "Unknown"}
YESNO = {"yes": "Yes", "no": "No", "unknown": "Unknown"}
OUTCOME = {"works": "Worked", "fails": "Did Not Work", "mixed": "Mixed", "untested": "Untested"}
VERDICT = {"extracellular": "Outside the Cell", "intracellular": "Inside the Cell", "mixed": "Spans Inside and Outside",
           "removed": "Removed From Mature Protein", "unknown": "Unknown"}
REGION = {"extracellular": "Outside the Cell", "transmembrane": "Membrane-Spanning", "cytoplasmic": "Inside the Cell",
          "lumenal": "Lumenal", "other": "Other"}
EVIDENCE_BASIS = {"search_verified": "Retrieved Source", "prior_knowledge": "General Knowledge", "hypothesis": "Hypothesis"}
EXPECTED = {"high": "High", "low": "Low", "absent": "Absent", "unknown": "Unknown"}
POLARITY = {"apical": "Apical", "basolateral": "Basolateral", "not_polarized": "Not Polarized", "unknown": "Unknown"}
APPLICATION = {"live_flow": "Live-Cell Flow Cytometry", "facs_sort": "Live-Cell Sort", "macs": "MACS Sort",
               "live_icc": "Live-Cell Imaging", "fixed_nonperm_flow": "Fixed Flow (Not Permeabilized)",
               "fixed_nonperm_icc": "Fixed Imaging (Not Permeabilized)", "perm_flow": "Permeabilized Flow",
               "perm_icc": "Permeabilized Imaging", "western_blot": "Western Blot", "ihc": "IHC"}
ROLE = {"CPEC_identity": "Keep: Marks CPECs", "pluripotent_contaminant": "Remove: Leftover Stem Cells",
        "neural_progenitor_contaminant": "Remove: Neural Progenitors",
        "endothelial_contaminant": "Remove: Endothelial Cells",
        "mesenchymal_contaminant": "Remove: Mesenchymal Cells", "pan_epithelial": "Check: All Epithelium"}
STAGE = {"candidate": "Candidate", "ordered": "Ordered", "received": "Received", "titrating": "Titrating",
         "testing": "Testing on Cells", "passed": "Passed", "failed": "Failed", "retired": "Retired"}
