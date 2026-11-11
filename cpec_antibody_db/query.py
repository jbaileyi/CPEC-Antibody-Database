"""Query helper for the CPEC antibody database.

Usage:
  python query.py gaps                 CPEC identity targets lacking a live-cell validated antibody
  python query.py live-ready           Antibodies with a successful test on non-permeabilized cells
  python query.py coverage             Antibody and validation counts per target
  python query.py panel [NAME]         Readiness of each panel member
  python query.py target SYMBOL        Full record for one target
  python query.py check                Integrity and consistency checks
  python query.py export OUTDIR        Write every table to CSV
  python query.py sql "SELECT ..."     Run a read-only query
"""
import csv
import sqlite3
import sys
from pathlib import Path

DB = Path(__file__).parent / "cpec_antibodies.sqlite"


def connect():
    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA foreign_keys = ON")
    return con


def show(rows):
    if not rows:
        print("(no rows)")
        return
    cols = rows[0].keys()
    data = [[("" if r[c] is None else str(r[c])) for c in cols] for r in rows]
    widths = [max(len(c), *(len(d[i]) for d in data)) for i, c in enumerate(cols)]
    widths = [min(w, 60) for w in widths]
    line = "  ".join(c.ljust(widths[i]) for i, c in enumerate(cols))
    print(line)
    print("  ".join("-" * w for w in widths))
    for d in data:
        print("  ".join(v[: widths[i]].ljust(widths[i]) for i, v in enumerate(d)))


def check(con):
    problems = []
    for row in con.execute("PRAGMA foreign_key_check"):
        problems.append(f"Foreign key violation: {tuple(row)}")
    rows = con.execute(
        """SELECT t.symbol, a.clone FROM validation v
           JOIN antibody a ON a.antibody_id = v.antibody_id
           JOIN target t ON t.target_id = a.target_id
           WHERE v.outcome = 'works'
             AND v.application IN ('live_flow','facs_sort','macs','live_icc')
             AND t.live_cell_accessible = 'no'"""
    ).fetchall()
    for r in rows:
        problems.append(f"{r['symbol']} {r['clone']}: live-cell success recorded for a target marked not accessible")
    rows = con.execute(
        """SELECT t.symbol FROM target t
           WHERE t.target_class IN ('single_pass_tm','multi_pass_tm','gpi_anchored')
             AND t.uniprot_acc IS NOT NULL AND t.uniprot_verified = 0"""
    ).fetchall()
    if rows:
        problems.append(
            "UniProt accession not yet verified (click Update in the app) for: " + ", ".join(r["symbol"] for r in rows)
        )
    rows = con.execute(
        """SELECT t.symbol FROM target t
           WHERE t.target_class IN ('single_pass_tm','multi_pass_tm')
             AND NOT EXISTS (SELECT 1 FROM target_region r
                              WHERE r.target_id = t.target_id AND r.region_type = 'extracellular')"""
    ).fetchall()
    if rows:
        problems.append("No extracellular region coordinates yet for: " + ", ".join(r["symbol"] for r in rows))
    rows = con.execute("SELECT COUNT(*) AS n FROM antibody WHERE curation_status = 'seed_unverified'").fetchone()
    if rows["n"]:
        problems.append(f"{rows['n']} antibody records are seed_unverified (no vendor or catalog number confirmed)")
    if problems:
        print("Findings:")
        for p in problems:
            print(" -", p)
    else:
        print("No problems found.")


def export(con, outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    tables = [r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type IN ('table','view') AND name NOT LIKE 'sqlite_%'")]
    for t in tables:
        rows = con.execute(f"SELECT * FROM {t}").fetchall()
        with open(out / f"{t}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            if rows:
                w.writerow(rows[0].keys())
                w.writerows([tuple(r) for r in rows])
    print(f"Wrote {len(tables)} files to {out}")


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return
    cmd = argv[1]
    con = connect()
    if cmd == "gaps":
        show(con.execute("SELECT * FROM v_gaps").fetchall())
    elif cmd == "live-ready":
        show(con.execute("SELECT * FROM v_live_ready").fetchall())
    elif cmd == "coverage":
        show(con.execute("SELECT * FROM v_coverage ORDER BY symbol").fetchall())
    elif cmd == "panel":
        if len(argv) > 2:
            show(con.execute("SELECT * FROM v_panel_readiness WHERE panel = ?", (argv[2],)).fetchall())
        else:
            show(con.execute("SELECT * FROM v_panel_readiness").fetchall())
    elif cmd == "target" and len(argv) > 2:
        sym = argv[2]
        for title, q in [
            ("Target", "SELECT * FROM target WHERE symbol = ? COLLATE NOCASE"),
            ("CPEC context", "SELECT c.* FROM cpec_context c JOIN target t USING(target_id) WHERE t.symbol = ? COLLATE NOCASE"),
            ("Regions", "SELECT r.* FROM target_region r JOIN target t USING(target_id) WHERE t.symbol = ? COLLATE NOCASE"),
            ("Antibodies", "SELECT a.* FROM antibody a JOIN target t USING(target_id) WHERE t.symbol = ? COLLATE NOCASE"),
            ("Validation", """SELECT v.* FROM validation v JOIN antibody a USING(antibody_id)
                              JOIN target t USING(target_id) WHERE t.symbol = ? COLLATE NOCASE"""),
        ]:
            print(f"\n== {title}")
            show(con.execute(q, (sym,)).fetchall())
    elif cmd == "check":
        check(con)
    elif cmd == "export" and len(argv) > 2:
        export(con, argv[2])
    elif cmd == "sql" and len(argv) > 2:
        if not argv[2].strip().lower().startswith("select"):
            raise SystemExit("Only SELECT statements are allowed here.")
        show(con.execute(argv[2]).fetchall())
    else:
        print(__doc__)


if __name__ == "__main__":
    main(sys.argv)
