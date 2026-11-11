"""Database access, first-run setup, and small migrations for older database files."""
import sqlite3
from pathlib import Path

HERE = Path(__file__).parent
DB = HERE / "cpec_antibodies.sqlite"

# Starting knowledge about how selective each seed marker is across tissues.
SPECIFICITY_SEED = {
    "ATP1A1": "ubiquitous", "ATP1B1": "ubiquitous",
    "EPCAM": "shared", "CD99": "shared", "FOLR1": "shared", "PECAM1": "shared", "NCAM1": "shared",
}

MIGRATION_SQL = """
CREATE TABLE IF NOT EXISTS literature_hit (
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
CREATE TABLE IF NOT EXISTS setting (
  key   TEXT PRIMARY KEY,
  value TEXT
);
CREATE TABLE IF NOT EXISTS sort_experiment (
  experiment_id INTEGER PRIMARY KEY,
  run_date TEXT, keep_markers TEXT, remove_markers TEXT, cells TEXT,
  input_cells TEXT, purity_before REAL, purity_after REAL, yield_pct REAL, viability_pct REAL,
  method TEXT, note TEXT
);
CREATE TABLE IF NOT EXISTS registry_summary (
  target_id INTEGER PRIMARY KEY REFERENCES target(target_id) ON DELETE CASCADE,
  n_records INTEGER, n_flow INTEGER, n_clonal INTEGER, n_conjugated INTEGER, top_vendor TEXT, fetched_on TEXT
);
CREATE TABLE IF NOT EXISTS target_info (
  target_id       INTEGER PRIMARY KEY REFERENCES target(target_id) ON DELETE CASCADE,
  function_text   TEXT, location_text TEXT, tissue_text TEXT, ptm_text TEXT,
  n_glycosylation INTEGER, n_disulfide INTEGER, n_isoforms INTEGER,
  shedding_hint   INTEGER NOT NULL DEFAULT 0,
  hpa_ensembl     TEXT, hpa_cp_ntpm TEXT, hpa_location TEXT, fetched_on TEXT
);
"""

INFO_COLUMNS = [("aliases", "TEXT"), ("mass_kda", "REAL"), ("n_pdb", "INTEGER"), ("length_aa", "INTEGER"),
                ("n_tm", "INTEGER"), ("has_signal", "INTEGER")]
ANTIBODY_COLUMNS = [("stage", "TEXT NOT NULL DEFAULT 'candidate'"), ("ordered_on", "TEXT"), ("received_on", "TEXT"),
                    ("lot", "TEXT"), ("price", "TEXT"), ("conjugate", "TEXT"), ("azide_free", "TEXT NOT NULL DEFAULT 'unknown'"),
                    ("flow_claimed", "TEXT NOT NULL DEFAULT 'unknown'"), ("applications_text", "TEXT"),
                    ("datasheet_url", "TEXT"), ("citations", "INTEGER"), ("source", "TEXT"),
                    ("shortlisted", "INTEGER NOT NULL DEFAULT 0")]


def connect():
    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA foreign_keys = ON")
    return con


def migrate(con):
    """Idempotent. Brings databases created by earlier versions up to date."""
    con.executescript(MIGRATION_SQL)
    cols = [r[1] for r in con.execute("PRAGMA table_info(target)")]
    if "specificity" not in cols:
        con.execute("ALTER TABLE target ADD COLUMN specificity TEXT NOT NULL DEFAULT 'unknown'")
    have = [r[1] for r in con.execute("PRAGMA table_info(antibody)")]
    for name, decl in ANTIBODY_COLUMNS:
        if name not in have:
            con.execute(f"ALTER TABLE antibody ADD COLUMN {name} {decl}")
    have_i = [r[1] for r in con.execute("PRAGMA table_info(target_info)")]
    for name, decl in INFO_COLUMNS:
        if name not in have_i:
            con.execute(f"ALTER TABLE target_info ADD COLUMN {name} {decl}")
    con.executescript("""CREATE TABLE IF NOT EXISTS titration (
      titration_id INTEGER PRIMARY KEY, antibody_id INTEGER NOT NULL REFERENCES antibody(antibody_id) ON DELETE CASCADE,
      run_date TEXT, dilution TEXT, mfi_pos REAL, mfi_neg REAL, sd_neg REAL, note TEXT);""")
    for sym, val in SPECIFICITY_SEED.items():
        con.execute("UPDATE target SET specificity = ? WHERE symbol = ? AND specificity = 'unknown'", (val, sym))
    con.commit()
    import catalog
    catalog.load_catalog(con)
    con.commit()


def ensure_db():
    if not DB.exists():
        con = sqlite3.connect(DB)
        con.executescript((HERE / "schema.sql").read_text())
        con.executescript((HERE / "seed.sql").read_text())
        con.commit()
        con.close()
    con = connect()
    migrate(con)
    con.close()


def get_setting(con, key, default=""):
    row = con.execute("SELECT value FROM setting WHERE key = ?", (key,)).fetchone()
    return row["value"] if row and row["value"] is not None else default


def set_setting(con, key, value):
    con.execute("INSERT INTO setting (key, value) VALUES (?, ?) "
                "ON CONFLICT(key) DO UPDATE SET value = excluded.value", (key, value))
    con.commit()
