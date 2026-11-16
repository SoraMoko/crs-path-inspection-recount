#!/usr/bin/env python3
"""Recount the values of Tables 3-7 from the anonymized observations in data/.

Only the standard library is used and no request is sent.
    python recount.py            print the tables and compare them with expected.json
"""
import csv
import json
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLS = (1, 2, 3, 4)
# G31 plus one change each; see README.md
VARIANTS = ("G31P", "G31_len4", "G31_cs_932237", "G31_cs_932237_932238", "G31_cs_31rules")


def read(name):
    with (HERE / "data" / name).open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def wilson(k, n):
    z = 1.959963984540054
    den = 1 + z * z / n
    centre = (k / n + z * z / (2 * n)) / den
    half = z * math.sqrt(k / n * (1 - k / n) / n + z * z / (4 * n * n)) / den
    return 100 * max(0.0, centre - half), 100 * min(1.0, centre + half)


def block_bootstrap(rows, hit, seed, reps=10000):
    """Percentile interval (%) when index blocks are resampled with replacement.

    `rows` are the paths that passed G0, in file order. The order of the blocks
    and the seed are kept as in the original analysis so that the published
    interval is obtained exactly.
    """
    blocks = {}
    for r in rows:
        b = blocks.setdefault(r["cluster"], [0, 0])
        b[0] += 1
        b[1] += hit(r)
    blocks = list(blocks.values())
    rng = random.Random(seed)
    rates = []
    for _ in range(reps):
        draw = [blocks[rng.randrange(len(blocks))] for _ in blocks]
        rates.append(100 * sum(b[1] for b in draw) / sum(b[0] for b in draw))
    rates.sort()
    return rates[int(reps * 0.025)], rates[int(reps * 0.975) - 1]


def r3(x):
    return round(x, 3)


# ---------------------------------------------------------------- Table 3, Section 3.1, Section 4.3.1
def main_sample():
    rows = read("main_sample.csv")
    by_pl = defaultdict(list)
    for r in rows:
        by_pl[int(r["pl"])].append(r)

    counts, old = {}, {}
    for pl in PLS:
        passed = [r for r in by_pl[pl] if r["G0"] == "200"]
        counts[str(pl)] = {
            "measured": len(by_pl[pl]),
            "G0_pass": len(passed),
            "G9_added": sum(r["G9"] == "403" for r in passed),
            "G31_added": sum(r["G31"] == "403" for r in passed),
        }
        old[str(pl)] = {
            "blocked_by_4_29_0": sum(r["v4_29_0"] == "403" for r in by_pl[pl]),
            "blocked_by_G0": sum(r["G0"] == "403" for r in by_pl[pl]),
            "decision_differs": sum(r["v4_29_0"] != r["G0"] for r in by_pl[pl]),
        }

    passed = [r for r in by_pl[3] if r["G0"] == "200"]
    n = len(passed)
    tests = {
        "G9_added": (lambda r: r["G9"] == "403", "zerobase-3-G9"),
        "G31_added": (lambda r: r["G31"] == "403", "zerobase-3-G31"),
        "G31_only": (lambda r: r["G31"] == "403" and r["G9"] == "200", "zerobase-3-only"),
        "G31P_added": (lambda r: r["G31P"] == "403", "zerobase-3-G31P"),
    }
    rates = {}
    for name, (hit, seed) in tests.items():
        k = sum(hit(r) for r in passed)
        lo, hi = wilson(k, n)
        blo, bhi = block_bootstrap(passed, hit, seed)
        rates[name] = {"count": k, "rate": r3(100 * k / n), "wilson": [r3(lo), r3(hi)], "bootstrap": [r3(blo), r3(bhi)]}

    pr = {}
    for pl in (3, 4):
        passed = [r for r in by_pl[pl] if r["G0"] == "200"]
        only31 = [r for r in passed if r["G31"] == "403" and r["G9"] == "200"]
        onlyp = [r for r in passed if r["G31P"] == "403" and r["G9"] == "200"]
        pr[str(pl)] = {
            "G31_only": len(only31),
            "G31P_only": len(onlyp),
            "G31_only_932237_alone": sum(r["path_rules_G31"] == "932237" for r in only31),
            "G31P_only_932237_alone": sum(r["path_rules_G31P"] == "932237" for r in onlyp),
        }
    return {"table3a": counts, "table3b": rates, "main_4_29_0_vs_G0": old, "main_G31P": pr}


# ---------------------------------------------------------------- Table 4
def removal_checks():
    rows = read("removal_checks.csv")
    out = {}
    for sample in ("main", "supplementary"):
        s = [r for r in rows if r["sample"] == sample]
        only = [r for r in s if r["G9"] == "200"]
        out[sample] = {
            "G31_added": len(s),
            "also_G9": len(s) - len(only),
            "G31_only": len(only),
            "pass_without_932237": sum(r["G31_minus_932237"] == "200" for r in only),
            "pass_without_932238": sum(r["G31_minus_932238"] == "200" for r in only),
            "pass_only_without_both": sum(r["G31_minus_932237"] == "403" and r["G31_minus_both"] == "200" for r in only),
            "still_blocked": sum(r["G31_minus_both"] == "403" for r in only),
        }
    out["responses"] = 5 * len(rows)
    return {"table4": out}


# ---------------------------------------------------------------- Section 4.1.5
def sites():
    rows = read("sites.csv")
    blocked = [r for r in rows if r["G31"] == "403"]
    return {"sites": {
        "sites": len(rows),
        "path_strings": len({r["path_string"] for r in rows}),
        "blocked_path_strings": len({r["path_string"] for r in blocked}),
        "blocked_by_G0_G9_4_29_0": [sum(r[g] == "403" for r in rows) for g in ("G0", "G9", "v4_29_0")],
        "pass_without_932237": sum(r["G31_minus_932237"] == "200" for r in blocked),
        "pass_without_932238": sum(r["G31_minus_932238"] == "200" for r in blocked),
        "pass_only_without_both": sum(r["G31_minus_932237"] == "403" and r["G31_minus_both"] == "200" for r in blocked),
        "blocked_by_G31P_minus_932238": sum(r["G31P_minus_932238"] == "403" for r in rows),
    }}


# ---------------------------------------------------------------- Table 5, Sections 3.1, 3.2.3, 4.2.1, 4.3.2
def full_group():
    rows = read("full_group.csv")
    cand = Counter(r["status"] for r in read("full_group_candidates.csv"))
    by_pl = defaultdict(list)
    for r in rows:
        by_pl[int(r["pl"])].append(r)
    t = {}
    for pl in PLS:
        s = by_pl[pl]
        only = [r for r in s if r["G9"] == "200" and r["G31"] == "403"]
        lo, hi = wilson(len(only), len(s))
        t[str(pl)] = {
            "inputs": len(s),
            "G0_blocked": sum(r["G0"] == "403" for r in s),
            "G9_blocked": sum(r["G9"] == "403" for r in s),
            "G31_blocked": sum(r["G31"] == "403" for r in s),
            "G31_only": len(only),
            "G31_only_pct": round(100 * len(only) / len(s), 1),
            "G31_only_wilson": [round(lo, 1), round(hi, 1)],
            "not_nested": sum((r["G0"] == "403" and r["G9"] == "200") or (r["G9"] == "403" and r["G31"] == "200") for r in s),
            "G31_only_by_path_score": sum(int(r["path_score_932_G31"]) >= 5 for r in only),
            "G31_only_with_932237": sum("932237" in r["path_rules_G31"].split() for r in only),
            "G31_only_932140_alone": sum(r["path_rules_G31"] == "932140" for r in only),
            "decision_differs_4_29_0_vs_G0": sum(r["v4_29_0"] != r["G0"] for r in s),
            "decision_differs_G31P_vs_G31": sum(r["G31P"] != r["G31"] for r in s),
            "G31P_blocked": sum(r["G31P"] == "403" for r in s),
        }
    return {"table5": t, "full_group_candidates": {"positive_tests": sum(cand.values()), "used": cand["used"],
                                                   "no_value": cand["no_value"], "too_long_or_empty": cand["too_long_or_empty"]}}


# ---------------------------------------------------------------- Section 4.2.1 (selected group)
def selected_group():
    rows = read("regression_inputs.csv")
    obs = {(r["id"], int(r["pl"]), r["config"]): r["observed"] for r in rows}
    ids = sorted({k[0] for k in obs})
    return {"selected_group": {
        "inputs": len(ids),
        "G31_only": [sum(obs[i, pl, "G9"] == "200" and obs[i, pl, "G31"] == "403" for i in ids) for pl in PLS],
    }}


# ---------------------------------------------------------------- Table 6
def detailed_inputs():
    rows = read("detailed_inputs.csv")
    http = defaultdict(set)
    for r in rows:
        http[r["id"], int(r["pl"]), r["config"], r["stage"]].add(r["http"])
    assert all(len(v) == 1 for v in http.values()), "replicates differ"
    http = {k: v.pop() for k, v in http.items()}
    source = {r["id"]: r["source_rule"] for r in rows}
    removed = {r["id"]: r["removed_path_check"] for r in rows if r["removed_path_check"]}

    def g31_only(i, stage):
        return [pl for pl in PLS
                if http[i, pl, "G31", stage] == "403" and http[i, pl, "G9", stage] == "200" and http[i, pl, "G0", stage] == "200"]

    out = {}
    for i in sorted(source):
        out[i] = {
            "source_rule": source[i],
            "G31_only_before_exclusion": g31_only(i, "before_exclusion"),
            "G31_only_after_exclusion": g31_only(i, "after_exclusion"),
            "needed_path_check_PL1": removed[i],
            "PL1_blocked_removed_restored": [http[i, 1, "G31", s] for s in ("after_exclusion", "path_check_removed", "path_check_restored")],
        }
    return {"table6": out, "table6_responses": len(rows)}


# ---------------------------------------------------------------- Table 7
def comparison():
    main = defaultdict(list)
    for r in read("main_sample.csv"):
        if r["G0"] == "200":
            main[int(r["pl"])].append(r)
    site = read("sites.csv")
    full = defaultdict(list)
    for r in read("full_group.csv"):
        full[int(r["pl"])].append(r)

    out = {}
    for g in ("G31",) + VARIANTS + ("G9",):
        out[g] = {
            "main_added_PL3": sum(r[g] == "403" for r in main[3]),
            "main_added_PL4": sum(r[g] == "403" for r in main[4]),
            "sites_blocked": sum(r[g] == "403" for r in site),
            "full_group_turned_to_pass": [sum(r["G31"] == "403" and r[g] == "200" for r in full[pl]) for pl in PLS],
            "newly_blocked": sum(r[g] == "403" and r["G31"] == "200" for pl in (3, 4) for r in main[pl])
                             + sum(r[g] == "403" and r["G31"] == "200" for pl in PLS for r in full[pl]),
        }
    out["G31"]["main_added_PL2"] = sum(r["G31"] == "403" for r in main[2])
    out["G31_cs_31rules"]["main_added_PL2"] = sum(r["G31_cs_31rules"] == "403" for r in main[2])
    lost = defaultdict(Counter)
    for g in ("G31_cs_932237", "G31_len4"):
        for r in full[3]:
            if r["G31"] == "403" and r[g] == "200":
                lost[g][r["source_rule"]] += 1
    return {"table7": out, "table7_notes": {
        "G31_blocked_PL3": sum(r["G31"] == "403" for r in full[3]),
        "kept_with_cs_932237_932238": sum(r["G31"] == "403" and r["G31_cs_932237_932238"] == "403" for r in full[3]),
        "turned_to_pass_by_source_rule": {g: dict(sorted(c.items())) for g, c in lost.items()},
    }}


def main():
    result = {}
    for part in (main_sample, removal_checks, sites, full_group, selected_group, detailed_inputs, comparison):
        result.update(part())
    print(json.dumps(result, indent=2))

    expected = json.loads((HERE / "expected.json").read_text(encoding="utf-8"))
    diff = [k for k in set(expected) | set(result) if expected.get(k) != result.get(k)]
    if diff:
        print("MISMATCH with expected.json:", ", ".join(sorted(diff)), file=sys.stderr)
        return 1
    print("all values match expected.json", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
