"""Panel designer checks: fluorochrome clashes, missing secondaries, azide, and gate logic."""
import re

CHANNELS = [  # (pattern, channel, laser)
    (r"fitc|alexa\s*fluor\s*488|af488|gfp", "FITC / 488", "blue"), (r"pe[-/ ]?cy7|pecy7", "PE-Cy7", "yellow-green"),
    (r"pe[-/ ]?cy5|pecy5", "PE-Cy5", "yellow-green"), (r"pe[-/ ]?dazzle|pe[-/ ]?texas|pe[-/ ]?cf594", "PE-Dazzle594", "yellow-green"),
    (r"\bpe\b|phycoerythrin", "PE", "yellow-green"), (r"percp", "PerCP", "blue"), (r"apc[-/ ]?cy7|apc[-/ ]?fire750", "APC-Cy7", "red"),
    (r"apc[-/ ]?r700|alexa\s*fluor\s*700|af700", "APC-R700 / 700", "red"), (r"alexa\s*fluor\s*647|af647|\bapc\b|allophycocyanin", "APC / 647", "red"),
    (r"bv421|brilliant violet 421|pacific blue", "BV421", "violet"), (r"bv510|brilliant violet 510", "BV510", "violet"),
    (r"bv605|brilliant violet 605", "BV605", "violet"), (r"bv650|brilliant violet 650", "BV650", "violet"),
    (r"bv711|brilliant violet 711", "BV711", "violet"), (r"bv786|brilliant violet 786", "BV786", "violet"),
    (r"alexa\s*fluor\s*594|af594", "AF594", "yellow-green"), (r"alexa\s*fluor\s*555|af555", "AF555", "yellow-green"),
    (r"biotin", "Biotin (Needs Streptavidin)", "none"),
]
# Channels that overlap strongly on common analyzers.
OVERLAP = [({"PE", "PE-Dazzle594"}, "PE and PE-Dazzle594 overlap heavily"), ({"APC / 647", "APC-R700 / 700"}, "APC and APC-R700 overlap"),
           ({"FITC / 488", "PerCP"}, "FITC and PerCP need compensation"), ({"PE-Cy5", "PerCP"}, "PE-Cy5 and PerCP overlap"),
           ({"PE", "AF555"}, "PE and AF555 are the same channel")]


def channel_of(conjugate):
    c = (conjugate or "").strip()
    if not c:
        return "Unconjugated", "none"
    for pat, name, laser in CHANNELS:
        if re.search(pat, c, re.I):
            return name, laser
    return c, "unknown"


def check_panel(members):
    """members: list of dicts with symbol, role ('keep' or 'remove'), conjugate, azide_free, verdict.
    Returns (rows, warnings)."""
    rows, warns, seen = [], [], {}
    for m in members:
        ch, laser = channel_of(m.get("conjugate"))
        rows.append({"Marker": m["symbol"], "Gate": "Keep" if m["role"] == "keep" else "Remove",
                     "Channel": ch, "Laser": laser, "Product": m.get("label", "")})
        seen.setdefault(ch, []).append(m["symbol"])
        if ch == "Unconjugated":
            warns.append(f"{m['symbol']}: unconjugated. It needs a labeled secondary antibody, which can bind other primaries from the same host.")
        if "Biotin" in ch:
            warns.append(f"{m['symbol']}: biotin conjugate needs a streptavidin step.")
        if m.get("azide_free") == "no":
            warns.append(f"{m['symbol']}: contains azide. Use an azide-free or LEAF version for cells you will culture.")
        if m.get("verdict") in ("intracellular", "removed"):
            warns.append(f"{m['symbol']}: immunogen is not outside the cell. It will not bind live cells.")
    for ch, syms in seen.items():
        if ch not in ("Unconjugated", "unknown") and len(syms) > 1:
            warns.append(f"Channel clash: {', '.join(syms)} all use {ch}. Give each marker its own channel.")
    chans = set(seen)
    for pair, msg in OVERLAP:
        if pair <= chans:
            warns.append(f"Spillover: {msg}.")
    if not any(m["role"] == "keep" for m in members):
        warns.append("No Keep marker chosen. A sort needs at least one marker that identifies CPECs.")
    if not any(m["role"] == "remove" for m in members):
        warns.append("No Remove marker chosen. Without one, leftover cells of other types stay in the sorted population.")
    warns.append("Always add a viability dye (for example DAPI or 7-AAD) in a channel not used above, plus FMO controls.")
    return rows, warns
