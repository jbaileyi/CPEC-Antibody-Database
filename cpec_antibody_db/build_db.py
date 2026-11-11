"""Build the CPEC extracellular-epitope antibody database from schema.sql and seed.sql.

Safe by default: if cpec_antibodies.sqlite already exists, this refuses to overwrite it, because the
file holds any records you added. Use --force to replace it. A timestamped backup is made first.
"""
import shutil
import sqlite3
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).parent
DB = HERE / "cpec_antibodies.sqlite"


def build(db_path=DB, with_seed=True, force=False):
    if db_path.exists():
        if not force:
            raise SystemExit(f"{db_path.name} already exists and may contain your records. "
                             "Re-run with --force to replace it (a backup is saved first).")
        backup = db_path.with_name(f"{db_path.stem}.backup-{datetime.now():%Y%m%d-%H%M%S}.sqlite")
        shutil.copy2(db_path, backup)
        print(f"Backup saved to {backup.name}")
        db_path.unlink()
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA foreign_keys = ON")
    con.executescript((HERE / "schema.sql").read_text())
    if with_seed:
        con.executescript((HERE / "seed.sql").read_text())
    con.commit()
    bad = con.execute("PRAGMA foreign_key_check").fetchall()
    con.close()
    if bad:
        raise SystemExit(f"Foreign key violations: {bad}")
    return db_path


if __name__ == "__main__":
    path = build(with_seed="--no-seed" not in sys.argv, force="--force" in sys.argv)
    print(f"Built {path}")
