-- Extracellular-epitope antibody database for stem-cell-derived choroid plexus epithelial cells (CPECs)
-- Design rule: a record only counts as evidence of extracellular accessibility if it was obtained
-- on non-permeabilized cells. The validation table enforces this with a CHECK constraint.

PRAGMA foreign_keys = ON;

CREATE TABLE reference (
  ref_id      INTEGER PRIMARY KEY,
  short_name  TEXT NOT NULL UNIQUE,
  citation    TEXT NOT NULL,
  doi         TEXT,
  pmid        TEXT,
  url         TEXT,
  note        TEXT
);

CREATE TABLE target (
  target_id            INTEGER PRIMARY KEY,
  symbol               TEXT NOT NULL,
  protein_name         TEXT,
  species              TEXT NOT NULL DEFAULT 'human' CHECK (species IN ('human','mouse','other')),
  uniprot_acc          TEXT,
  uniprot_verified     INTEGER NOT NULL DEFAULT 0 CHECK (uniprot_verified IN (0,1)),
  target_class         TEXT NOT NULL CHECK (target_class IN
                         ('single_pass_tm','multi_pass_tm','gpi_anchored','glycan_epitope','glycoprotein_epitope')),
  live_cell_accessible TEXT NOT NULL DEFAULT 'unknown' CHECK (live_cell_accessible IN ('yes','partial','no','unknown')),
  accessibility_note   TEXT,
  specificity          TEXT NOT NULL DEFAULT 'unknown' CHECK (specificity IN ('selective','shared','ubiquitous','unknown')),
  UNIQUE (symbol, species)
);

-- Extracellular, transmembrane, and cytoplasmic segments. Populated by ingest_uniprot.py.
CREATE TABLE target_region (
  region_id    INTEGER PRIMARY KEY,
  target_id    INTEGER NOT NULL REFERENCES target(target_id) ON DELETE CASCADE,
  region_type  TEXT NOT NULL CHECK (region_type IN ('extracellular','transmembrane','cytoplasmic','lumenal','other')),
  start_aa     INTEGER,
  end_aa       INTEGER,
  source       TEXT NOT NULL,
  retrieved_on TEXT,
  note         TEXT
);

-- What role a target plays when sorting or characterizing derived CPECs.
CREATE TABLE cpec_context (
  context_id         INTEGER PRIMARY KEY,
  target_id          INTEGER NOT NULL REFERENCES target(target_id) ON DELETE CASCADE,
  context            TEXT NOT NULL CHECK (context IN
                       ('CPEC_identity','pluripotent_contaminant','neural_progenitor_contaminant',
                        'endothelial_contaminant','mesenchymal_contaminant','pan_epithelial')),
  expected_in_cpec   TEXT NOT NULL DEFAULT 'unknown' CHECK (expected_in_cpec IN ('high','low','absent','unknown')),
  polarity           TEXT NOT NULL DEFAULT 'unknown' CHECK (polarity IN ('apical','basolateral','not_polarized','unknown')),
  stage_note         TEXT,
  evidence_basis     TEXT NOT NULL CHECK (evidence_basis IN ('search_verified','prior_knowledge','hypothesis')),
  ref_id             INTEGER REFERENCES reference(ref_id),
  note               TEXT,
  UNIQUE (target_id, context)
);

CREATE TABLE epitope (
  epitope_id      INTEGER PRIMARY KEY,
  target_id       INTEGER NOT NULL REFERENCES target(target_id) ON DELETE CASCADE,
  epitope_type    TEXT NOT NULL CHECK (epitope_type IN ('linear','conformational','glycan','glycoprotein_dependent','unknown')),
  label           TEXT,
  region_type     TEXT,
  start_aa        INTEGER,
  end_aa          INTEGER,
  glycan_dependent INTEGER NOT NULL DEFAULT 0 CHECK (glycan_dependent IN (0,1)),
  mapping_method  TEXT,
  note            TEXT
);

CREATE TABLE antibody (
  antibody_id     INTEGER PRIMARY KEY,
  target_id       INTEGER NOT NULL REFERENCES target(target_id) ON DELETE CASCADE,
  epitope_id      INTEGER REFERENCES epitope(epitope_id),
  clone           TEXT,
  host            TEXT,
  isotype         TEXT,
  clonality       TEXT CHECK (clonality IN ('monoclonal','polyclonal','recombinant','unknown')),
  reactivity      TEXT,
  vendor          TEXT,
  catalog_no      TEXT,
  rrid            TEXT,
  immunogen       TEXT,
  conjugate       TEXT,
  azide_free      TEXT NOT NULL DEFAULT 'unknown',
  flow_claimed    TEXT NOT NULL DEFAULT 'unknown',
  applications_text TEXT,
  datasheet_url   TEXT,
  citations       INTEGER,
  source          TEXT,
  shortlisted     INTEGER NOT NULL DEFAULT 0,
  curation_status TEXT NOT NULL DEFAULT 'needs_curation' CHECK (curation_status IN
                    ('seed_unverified','needs_curation','vendor_verified','lab_verified')),
  note            TEXT
);

-- One row per application test. Any application that reads the extracellular epitope must be on
-- non-permeabilized cells, so the constraint below blocks inconsistent entries.
CREATE TABLE validation (
  validation_id        INTEGER PRIMARY KEY,
  antibody_id          INTEGER NOT NULL REFERENCES antibody(antibody_id) ON DELETE CASCADE,
  application          TEXT NOT NULL CHECK (application IN
                         ('live_flow','facs_sort','macs','live_icc','fixed_nonperm_flow','fixed_nonperm_icc',
                          'perm_flow','perm_icc','western_blot','ihc')),
  cell_context         TEXT,
  species              TEXT,
  context_matches_cpec INTEGER NOT NULL DEFAULT 0 CHECK (context_matches_cpec IN (0,1)),
  nonpermeabilized     INTEGER NOT NULL DEFAULT 0 CHECK (nonpermeabilized IN (0,1)),
  positive_control     TEXT,
  negative_control     TEXT,
  ko_or_kd_validated   INTEGER NOT NULL DEFAULT 0 CHECK (ko_or_kd_validated IN (0,1)),
  outcome              TEXT NOT NULL CHECK (outcome IN ('works','fails','mixed','untested')),
  ref_id               INTEGER REFERENCES reference(ref_id),
  note                 TEXT,
  CHECK (application NOT IN ('live_flow','facs_sort','macs','live_icc','fixed_nonperm_flow','fixed_nonperm_icc')
         OR nonpermeabilized = 1)
);

CREATE TABLE panel (
  panel_id    INTEGER PRIMARY KEY,
  name        TEXT NOT NULL UNIQUE,
  purpose     TEXT,
  stage       TEXT,
  status      TEXT NOT NULL DEFAULT 'draft_hypothesis' CHECK (status IN ('draft_hypothesis','tested','retired')),
  note        TEXT
);

-- Panels reference targets first so a panel can be designed before a validated antibody exists.
CREATE TABLE panel_member (
  panel_member_id INTEGER PRIMARY KEY,
  panel_id        INTEGER NOT NULL REFERENCES panel(panel_id) ON DELETE CASCADE,
  target_id       INTEGER NOT NULL REFERENCES target(target_id),
  antibody_id     INTEGER REFERENCES antibody(antibody_id),
  selection       TEXT NOT NULL CHECK (selection IN ('positive','negative')),
  fluorochrome    TEXT,
  note            TEXT
);

CREATE INDEX idx_antibody_target   ON antibody(target_id);
CREATE INDEX idx_validation_ab     ON validation(antibody_id);
CREATE INDEX idx_context_target    ON cpec_context(target_id);
CREATE INDEX idx_region_target     ON target_region(target_id);

-- Antibodies with a successful test on non-permeabilized cells.
CREATE VIEW v_live_ready AS
SELECT a.antibody_id, t.symbol, a.clone, a.host, a.isotype, a.curation_status,
       v.application, v.cell_context, v.context_matches_cpec, v.ko_or_kd_validated,
       r.short_name AS reference
FROM validation v
JOIN antibody a ON a.antibody_id = v.antibody_id
JOIN target   t ON t.target_id   = a.target_id
LEFT JOIN reference r ON r.ref_id = v.ref_id
WHERE v.outcome = 'works'
  AND v.application IN ('live_flow','facs_sort','macs','live_icc');

-- Coverage per target.
CREATE VIEW v_coverage AS
SELECT t.target_id, t.symbol, t.live_cell_accessible,
       (SELECT group_concat(c.context, ', ') FROM cpec_context c WHERE c.target_id = t.target_id) AS contexts,
       (SELECT COUNT(*) FROM antibody a WHERE a.target_id = t.target_id) AS n_antibodies,
       (SELECT COUNT(DISTINCT a.antibody_id) FROM antibody a
          JOIN validation v ON v.antibody_id = a.antibody_id
         WHERE a.target_id = t.target_id AND v.outcome = 'works'
           AND v.application IN ('live_flow','facs_sort','macs','live_icc')) AS n_live_validated,
       (SELECT COUNT(DISTINCT a.antibody_id) FROM antibody a
          JOIN validation v ON v.antibody_id = a.antibody_id
         WHERE a.target_id = t.target_id AND v.outcome = 'works'
           AND v.application IN ('live_flow','facs_sort','macs','live_icc')
           AND v.context_matches_cpec = 1) AS n_live_validated_in_cpec
FROM target t;

-- CPEC identity targets with no live-cell validated antibody in a CPEC context.
CREATE VIEW v_gaps AS
SELECT c.symbol, c.live_cell_accessible, c.n_antibodies, c.n_live_validated, c.n_live_validated_in_cpec
FROM v_coverage c
WHERE EXISTS (SELECT 1 FROM cpec_context x WHERE x.target_id = c.target_id AND x.context = 'CPEC_identity')
  AND c.n_live_validated_in_cpec = 0
ORDER BY c.n_antibodies ASC, c.symbol;

-- Readiness of each panel member.
CREATE VIEW v_panel_readiness AS
SELECT p.name AS panel, t.symbol, pm.selection, a.clone AS assigned_clone,
       (SELECT COUNT(*) FROM validation v
         WHERE v.antibody_id = pm.antibody_id AND v.outcome = 'works'
           AND v.application IN ('live_flow','facs_sort','macs','live_icc')) AS assigned_ab_live_validated,
       (SELECT COUNT(DISTINCT a2.antibody_id) FROM antibody a2
          JOIN validation v ON v.antibody_id = a2.antibody_id
         WHERE a2.target_id = t.target_id AND v.outcome = 'works'
           AND v.application IN ('live_flow','facs_sort','macs','live_icc')) AS any_ab_live_validated_on_target
FROM panel_member pm
JOIN panel  p ON p.panel_id  = pm.panel_id
JOIN target t ON t.target_id = pm.target_id
LEFT JOIN antibody a ON a.antibody_id = pm.antibody_id;

-- Papers found automatically from Europe PMC (leads to read, not verified evidence).
CREATE TABLE literature_hit (
  hit_id     INTEGER PRIMARY KEY,
  target_id  INTEGER NOT NULL REFERENCES target(target_id) ON DELETE CASCADE,
  source_id  TEXT NOT NULL,
  title      TEXT,
  year       TEXT,
  journal    TEXT,
  doi        TEXT,
  fetched_on TEXT,
  UNIQUE (target_id, source_id)
);

-- Per-user settings such as author name and default cell description.
CREATE TABLE setting (
  key   TEXT PRIMARY KEY,
  value TEXT
);

-- Facts about each marker, filled automatically from UniProt and the Human Protein Atlas.
CREATE TABLE target_info (
  target_id       INTEGER PRIMARY KEY REFERENCES target(target_id) ON DELETE CASCADE,
  function_text   TEXT,
  location_text   TEXT,
  tissue_text     TEXT,
  ptm_text        TEXT,
  n_glycosylation INTEGER,
  n_disulfide     INTEGER,
  n_isoforms      INTEGER,
  shedding_hint   INTEGER NOT NULL DEFAULT 0,
  hpa_ensembl     TEXT,
  hpa_cp_ntpm     TEXT,
  hpa_location    TEXT,
  fetched_on      TEXT
);
