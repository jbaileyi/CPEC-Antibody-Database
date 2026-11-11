"""Extra starting markers and clone names, loaded once into any database (safe to run again).

Everything here is general knowledge, not looked up in this build, and is labeled that way in the database:
roles use evidence basis 'prior_knowledge', clones are 'seed_unverified' (clone name only, no vendor or catalog
number confirmed). Topology, facts and RNA levels are filled in by 'Update From Public Databases'.
"""

CATALOG_VERSION = "2"

# symbol, protein name, class, role, expected in CPECs, polarity, specificity, note
MARKERS = [
    # Keep: markers reported or expected on choroid plexus epithelium
    ("LRP2", "Low-density lipoprotein receptor-related protein 2 (megalin)", "single_pass_tm", "CPEC_identity", "high", "apical", "shared",
     "Endocytic receptor reported in choroid plexus epithelium. Also made by kidney proximal tubule and other epithelia."),
    ("SLC2A1", "Solute carrier family 2, member 1 (GLUT1)", "multi_pass_tm", "CPEC_identity", "high", "basolateral", "ubiquitous",
     "Glucose transporter reported in choroid plexus. Made by many cells, so it needs a partner marker. One large extracellular loop."),
    ("TFRC", "Transferrin receptor protein 1 (CD71)", "single_pass_tm", "CPEC_identity", "high", "basolateral", "ubiquitous",
     "Iron uptake receptor, reported in choroid plexus. Also high on dividing cells, so it does not separate CPECs from progenitors."),
    ("SLC4A10", "Sodium-driven chloride bicarbonate exchanger (NCBE)", "multi_pass_tm", "CPEC_identity", "high", "basolateral", "shared",
     "Reported on the basolateral membrane of choroid plexus epithelium. Extracellular loops are short, so antibody choice is hard."),
    ("KCNJ13", "Inward rectifier potassium channel Kir7.1", "multi_pass_tm", "CPEC_identity", "high", "apical", "shared",
     "Reported as a choroid plexus epithelium channel. Two transmembrane helices with one extracellular loop."),
    ("CLDN2", "Claudin-2", "multi_pass_tm", "CPEC_identity", "high", "unknown", "shared",
     "Tight junction protein reported in choroid plexus. Extracellular loops may be hidden in mature junctions."),
    ("CLDN3", "Claudin-3", "multi_pass_tm", "CPEC_identity", "unknown", "unknown", "shared",
     "Tight junction protein of many epithelia. Possible epithelial identity marker."),
    ("SLC16A1", "Monocarboxylate transporter 1 (MCT1)", "multi_pass_tm", "CPEC_identity", "unknown", "unknown", "ubiquitous",
     "Widely made. Little extracellular sequence."),
    ("BSG", "Basigin (CD147)", "single_pass_tm", "CPEC_identity", "unknown", "unknown", "ubiquitous",
     "Chaperone for monocarboxylate transporters. Widely made, with a large glycosylated extracellular domain."),
    # Epithelial identity (all epithelium)
    ("CDH1", "Cadherin-1 (E-cadherin, CD324)", "single_pass_tm", "pan_epithelial", "high", "basolateral", "shared",
     "Epithelial adhesion protein. Confirms epithelial identity. Adhesion function can mean harsh detachment is needed before staining."),
    ("ITGA6", "Integrin alpha-6 (CD49f)", "single_pass_tm", "pan_epithelial", "unknown", "basolateral", "shared",
     "Epithelial and stem cell surface marker. Also found on leftover stem cells."),
    ("CLDN4", "Claudin-4", "multi_pass_tm", "pan_epithelial", "unknown", "unknown", "shared", "Tight junction protein of many epithelia."),
    # Remove: leftover pluripotent stem cells
    ("TNFRSF8", "Tumor necrosis factor receptor superfamily member 8 (CD30)", "single_pass_tm", "pluripotent_contaminant", "absent", "unknown", "shared",
     "Reported on human pluripotent stem cells. Also found on some lymphocytes."),
    ("ALPL", "Alkaline phosphatase, tissue-nonspecific (TNAP)", "gpi_anchored", "pluripotent_contaminant", "absent", "unknown", "shared",
     "GPI-anchored enzyme reported on pluripotent stem cells (the TRA-2-54 epitope). Also made in bone and other tissues."),
    ("CD9", "CD9 antigen", "multi_pass_tm", "pluripotent_contaminant", "unknown", "unknown", "ubiquitous",
     "Tetraspanin reported high on pluripotent stem cells but made by many cells, so a weak remover on its own. Small extracellular loops."),
    ("THY1", "Thy-1 membrane glycoprotein (CD90)", "gpi_anchored", "mesenchymal_contaminant", "absent", "unknown", "shared",
     "GPI-anchored marker of mesenchymal cells and fibroblasts. Also on neurons and some stem cells."),
    # Remove: neural progenitors
    ("FUT4", "Alpha-(1,3)-fucosyltransferase 4 (makes the SSEA-1 or CD15 epitope)", "single_pass_tm", "neural_progenitor_contaminant", "unknown", "unknown", "shared",
     "The CD15 (SSEA-1) glycan epitope is reported on neural stem cells and absent from naive human pluripotent cells. The antibody detects the glycan, not the enzyme."),
    ("CXCR4", "C-X-C chemokine receptor type 4 (CD184)", "multi_pass_tm", "neural_progenitor_contaminant", "unknown", "unknown", "shared",
     "Reported on neural progenitors and several other lineages. Extracellular N-terminus and loops."),
    ("NGFR", "Tumor necrosis factor receptor superfamily member 16 (CD271)", "single_pass_tm", "neural_progenitor_contaminant", "unknown", "unknown", "shared",
     "Reported on neural crest and mesenchymal progenitors."),
    ("L1CAM", "Neural cell adhesion molecule L1 (CD171)", "single_pass_tm", "neural_progenitor_contaminant", "unknown", "unknown", "shared",
     "Large extracellular domain. Reported on neurons and neural progenitors, and also on some epithelia."),
    ("SLC1A3", "Excitatory amino acid transporter 1 (GLAST)", "multi_pass_tm", "neural_progenitor_contaminant", "unknown", "unknown", "shared",
     "Radial glia and astrocyte marker. Extracellular epitopes are limited."),
    ("PDGFRA", "Platelet-derived growth factor receptor alpha (CD140a)", "single_pass_tm", "neural_progenitor_contaminant", "unknown", "unknown", "shared",
     "Reported on neural progenitors and mesenchymal cells. Large extracellular domain."),
    # Remove: endothelial
    ("CDH5", "Cadherin-5 (VE-cadherin, CD144)", "single_pass_tm", "endothelial_contaminant", "absent", "unknown", "selective",
     "Endothelial adhesion protein."),
    ("KDR", "Vascular endothelial growth factor receptor 2 (CD309)", "single_pass_tm", "endothelial_contaminant", "absent", "unknown", "shared",
     "Reported on endothelial and early mesoderm cells."),
    ("CD34", "Hematopoietic progenitor cell antigen CD34", "single_pass_tm", "endothelial_contaminant", "absent", "unknown", "shared",
     "Reported on endothelial and hematopoietic progenitors. Heavily glycosylated."),
    ("ENG", "Endoglin (CD105)", "single_pass_tm", "endothelial_contaminant", "absent", "unknown", "shared",
     "Reported on endothelial cells and mesenchymal cells."),
    # Remove: mesenchymal
    ("PDGFRB", "Platelet-derived growth factor receptor beta (CD140b)", "single_pass_tm", "mesenchymal_contaminant", "absent", "unknown", "shared",
     "Reported on pericytes and mesenchymal cells."),
    ("NT5E", "5'-nucleotidase (CD73)", "gpi_anchored", "mesenchymal_contaminant", "unknown", "unknown", "shared",
     "GPI-anchored enzyme reported on mesenchymal cells. Also found on some epithelia, so test it on CPECs."),
    ("CD44", "CD44 antigen", "single_pass_tm", "mesenchymal_contaminant", "unknown", "unknown", "ubiquitous",
     "Widely made, many splice forms. Weak remover on its own."),
]

# symbol, clone, host, isotype, conjugate-free note. Clone names only; vendor and catalog stay blank until confirmed.
CLONES = [
    ("EPCAM", "9C4", "Mouse", "IgG2b"), ("PECAM1", "WM59", "Mouse", "IgG1"), ("CD34", "581", "Mouse", "IgG1"),
    ("THY1", "5E10", "Mouse", "IgG1"), ("NT5E", "AD2", "Mouse", "IgG1"), ("ENG", "43A3", "Mouse", "IgG1"),
    ("CDH5", "BV9", "Mouse", "IgG1"), ("KDR", "7D4-6", "Mouse", "IgG1"), ("CXCR4", "12G5", "Mouse", "IgG2a"),
    ("NGFR", "ME20.4", "Mouse", "IgG1"), ("CD9", "HI9a", "Mouse", "IgG1"), ("TFRC", "CY1G4", "Mouse", "IgG2a"),
    ("BSG", "HIM6", "Mouse", "IgG1"), ("PDGFRB", "18A1", "Mouse", "IgG1"), ("ITGA6", "GoH3", "Rat", "IgG2a"),
    ("FUT4", "HI98", "Mouse", "IgM"), ("TNFRSF8", "BerH8", "Mouse", "IgG1"), ("ALPL", "W8B2", "Mouse", "IgG1"),
    ("CDH1", "67A4", "Mouse", "IgG1"), ("PROM1", "293C3", "Mouse", "IgG2b"), ("PDGFRA", "16A1", "Mouse", "IgG1"),
    ("L1CAM", "5G3", "Mouse", "IgG1"),
]


def load_catalog(con):
    """Adds the markers and clone names above if they are not already present. Never changes existing rows."""
    from db import get_setting, set_setting
    if get_setting(con, "catalog_version") == CATALOG_VERSION:
        return 0
    added = 0
    for sym, name, klass, role, exp, pol, spec, note in MARKERS:
        if con.execute("SELECT 1 FROM target WHERE symbol = ? AND species = 'human'", (sym,)).fetchone():
            tid = con.execute("SELECT target_id FROM target WHERE symbol = ? AND species = 'human'", (sym,)).fetchone()[0]
        else:
            con.execute("INSERT INTO target (symbol, protein_name, species, target_class, live_cell_accessible, "
                        "accessibility_note, specificity) VALUES (?,?,?,?,?,?,?)",
                        (sym, name, "human", klass, "unknown",
                         "Topology not loaded yet. Click Update From Public Databases.", spec))
            tid = con.execute("SELECT last_insert_rowid()").fetchone()[0]
            added += 1
        con.execute("INSERT OR IGNORE INTO cpec_context (target_id, context, expected_in_cpec, polarity, evidence_basis, note) "
                    "VALUES (?,?,?,?,?,?)",
                    (tid, role, exp, pol, "prior_knowledge", note + " General knowledge, not re-verified in this build."))
    for sym, clone, host, iso in CLONES:
        row = con.execute("SELECT target_id FROM target WHERE symbol = ? AND species = 'human'", (sym,)).fetchone()
        if not row:
            continue
        if con.execute("SELECT 1 FROM antibody WHERE target_id = ? AND clone = ?", (row[0], clone)).fetchone():
            continue
        con.execute("INSERT INTO antibody (target_id, clone, host, isotype, reactivity, curation_status, source, note) "
                    "VALUES (?,?,?,?,?,?,?,?)",
                    (row[0], clone, host, iso, "human", "seed_unverified", "Starting list",
                     "Clone name from general knowledge. Vendor, catalog number and datasheet not checked."))
    set_setting(con, "catalog_version", CATALOG_VERSION)
    return added
