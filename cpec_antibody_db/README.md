# CPEC antibody database

Finds and tracks antibodies that bind the outside of live cells, for sorting stem-cell-derived choroid plexus epithelial cells (CPECs).

## Start (no coding)

1. Install Python 3 from https://www.python.org/downloads/ (Windows: tick "Add Python to PATH").
2. Unzip this folder.
3. Double-click **run.bat** (Windows) or run **./run.sh** (Mac or Linux).
4. Your browser opens the database. The first launch installs what it needs, about a minute.

## Your routine

1. **Test next**: click *Update from public databases*, then open the top marker.
2. **Find antibodies**: search the Antibody Registry for that marker, read the datasheet checks, add products to the order list.
3. Run the staining plan on live, unpermeabilized cells.
4. Back on **Test next**, record *Worked*, *Did not work* or *Mixed*.

## What it does automatically

- **Topology and facts** (UniProt): which parts of the protein are outside the cell, function, location, tissue, glycosylation sites, disulfide bonds, isoforms, and whether the protein may be shed.
- **Choroid plexus RNA** (Human Protein Atlas).
- **Vendor lookup** (Antibody Registry): vendor, catalog number, RRID, clone, conjugate, isotype, listed applications and citation counts, plus search links for BioLegend, CiteAb and Biocompare.
- **Datasheet checks** for every product: immunogen outside the cell, flow listing, azide, conjugate, clonality, citations.
- **Ranking** of what to test next, a draft sort strategy, a staining plan with controls, and an order list you can download.
- **Marker discovery** (Human Protein Atlas): lists plasma membrane proteins with high choroid plexus RNA and adds the ones you pick, with topology downloaded. RNA is a lead, not proof of surface protein or selectivity.
- **Evidence labels**: every marker shows only what is recorded for it (your tests, live-cell sorts reported elsewhere, or the basis for its link to CPECs). Nothing is assumed.
- **Panel Designer**: checks fluorochrome clashes, unconjugated antibodies, azide and immunogen location for a Keep and Remove panel.
- **Antibody Tracker**: stage, order dates, lot, price, titration log with staining index.
- **Experiments**: log sort runs and see average purity gain, yield and viability.
- **Overview**: coverage, registry counts per marker, and gaps to fill.
- **Papers** (Europe PMC) that mention a marker with choroid plexus and sorting terms.

## Limits

- The Registry lists what vendors submit, so entries can be incomplete. Azide status is not in it and shows as unknown until you check the datasheet. Confirm catalog numbers and prices with the vendor.
- Accessibility estimates and scores are rules of thumb. Only your own test shows whether an antibody works on your cells.
- Nothing here publishes anything or creates a DOI.

## Your data

Everything is saved in `cpec_antibodies.sqlite`. Back it up. *My Records* has downloads (CSV zip and the database file) if you want to deposit the data somewhere yourself, for example Zenodo, which issues its own DOI. `build_db.py --force` starts over (it saves a backup first). `test_logic.py` checks the automatic parts without internet.

## Files

`app.py` interface, `autofill.py` UniProt, Protein Atlas, Europe PMC and Antibody Registry, `ranking.py` scoring and datasheet checks, `db.py` storage, `schema.sql` and `seed.sql` starting data, `query.py` command line queries.

## Seed data

47 markers, 33 antibody records (most are clone names only, `seed_unverified`), 8 references. Added markers and clones are general knowledge, not looked up. Roles of markers come from patents, preprints, a CIRM award summary and general knowledge, labeled by evidence level.
