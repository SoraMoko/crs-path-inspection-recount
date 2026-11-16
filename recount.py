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


# ---------------------------------------------------------------- Table 3, Table 7(a)
def main_sample():
    rows = read("main_sample.csv")
    by_pl = defaultdict(list)
    for r in rows:
        by_pl[int(r["pl"])].append(r)

    counts = {}
    for pl in PLS:
        passed = [r for r in by_pl[pl] if r["G0"] == "200"]
        counts[str(pl)] = {
            "measured": len(by_pl[pl]),
            "G0_pass": len(passed),
            "G9_added": sum(r["G9"] == "403" for r in passed),
            "G31_added": sum(r["G31"] == "403" for r in passed),
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
            "G31_added": sum(r["G31"] == "403" for r in passed),
            "G31P_added": sum(r["G31P"] == "403" for r in passed),
            "G31_only": len(only31),
            "G31P_only": len(onlyp),
            "G31_only_932237_alone": sum(r["path_rules_G31"] == "932237" for r in only31),
            "G31P_only_932237_alone": sum(r["path_rules_G31P"] == "932237" for r in onlyp),
            "newly_blocked_by_G31P": sum(r["G31P"] == "403" and r["G31"] == "200" for r in passed),
        }
    return {"table3a": counts, "table3b": rates, "table7a": pr}


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


# ---------------------------------------------------------------- Section 4.1.5, Table 7(b)
def sites():
    rows = read("sites.csv")
    blocked = [r for r in rows if r["G31"] == "403"]
    return {
        "sites": {
            "sites": len(rows),
            "path_strings": len({r["path_string"] for r in rows}),
            "blocked_path_strings": len({r["path_string"] for r in blocked}),
            "G0_blocked": sum(r["G0"] == "403" for r in rows),
            "G9_blocked": sum(r["G9"] == "403" for r in rows),
            "pass_without_932237": sum(r["G31_minus_932237"] == "200" for r in blocked),
            "pass_without_932238": sum(r["G31_minus_932238"] == "200" for r in blocked),
            "pass_only_without_both": sum(r["G31_minus_932237"] == "403" and r["G31_minus_both"] == "200" for r in blocked),
        },
        "table7b": {
            "G31": [len(blocked), len(rows) - len(blocked)],
            "G31P": [sum(r["G31P"] == "403" for r in rows), sum(r["G31P"] == "200" for r in rows)],
            "G31P_minus_932238": [sum(r["G31P_minus_932238"] == "403" for r in rows), sum(r["G31P_minus_932238"] == "200" for r in rows)],
        },
    }


# ---------------------------------------------------------------- Table 5
def regression_inputs():
    rows = read("regression_inputs.csv")
    out = {}
    for group in ("selected", "random"):
        obs, pred = {}, {}
        for r in rows:
            if r["group"] == group:
                obs[r["id"], int(r["pl"]), r["config"]] = r["observed"]
                pred[r["id"], int(r["pl"]), r["config"]] = r["predicted"]
        ids = sorted({k[0] for k in obs})
        t = {}
        for pl in PLS:
            t[str(pl)] = {
                "G9_blocked": sum(obs[i, pl, "G9"] == "403" for i in ids),
                "G31_blocked": sum(obs[i, pl, "G31"] == "403" for i in ids),
                "G31_only": sum(obs[i, pl, "G9"] == "200" and obs[i, pl, "G31"] == "403" for i in ids),
                "G31_only_predicted": sum(pred[i, pl, "G9"] == "200" and pred[i, pl, "G31"] == "403" for i in ids),
                "prediction_agrees": sum(obs[i, pl, g] == pred[i, pl, g] for i in ids for g in ("G9", "G31")),
            }
        agree = sum(v["prediction_agrees"] for v in t.values())
        wrong = [k for k in obs if obs[k] != pred[k]]
        out[group] = {
            "inputs": len(ids),
            "by_pl": t,
            "conditions": len(obs),
            "prediction_agrees": agree,
            "predicted_pass_observed_block": sum(pred[k] == "200" and obs[k] == "403" for k in wrong),
        }
    return {"table5": out}


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


# ---------------------------------------------------------------- Section 3.2.3
def random_group_candidates():
    rows = read("random_group_candidates.csv")
    c = Counter(r["status"] for r in rows)

    # A candidate is dropped when it is the same test as, or gives the same path as, a selected input.
    selected = [r for r in read("regression_inputs.csv") if r["group"] == "selected"]
    tests = {(r["source_rule"], r["test_id"]) for r in selected}
    keys = {r["path_key"] for r in selected}

    def status(r):
        if not r["path_key"]:
            return r["status"]
        if (r["source_rule"], r["test_id"]) in tests:
            return "same_test_as_selected"
        return "same_path_as_selected" if r["path_key"] in keys else "pool"

    pool = sorted((r for r in rows if r["status"] == "pool"), key=lambda r: (r["source_rule"], int(r["test_id"])))
    drawn = random.Random(20261109).sample(pool, 24)
    same_draw = [(r["source_rule"], r["test_id"]) for r in drawn] == \
        [(r["source_rule"], r["test_id"]) for r in sorted((r for r in rows if r["drawn_as"]), key=lambda r: r["drawn_as"])]
    return {"random_group_candidates": {
        "positive_tests": len(rows),
        "value_extracted": len(rows) - c["no_value"] - c["too_long_or_empty"],
        "overlap_with_selected": c["same_test_as_selected"] + c["same_path_as_selected"],
        "pool": len(pool),
        "exclusion_reproduced": all(status(r) == r["status"] for r in rows),
        "drawn": sum(bool(r["drawn_as"]) for r in rows),
        "draw_reproduced": same_draw,
    }}


def main():
    result = {}
    for part in (main_sample, removal_checks, sites, regression_inputs, detailed_inputs, random_group_candidates):
        result.update(part())
    print(json.dumps(result, indent=2))

    expected = json.loads((HERE / "expected.json").read_text(encoding="utf-8"))
    diff = [k for k in expected if expected[k] != result.get(k)]
    if diff:
        print("MISMATCH with expected.json:", ", ".join(diff), file=sys.stderr)
        return 1
    print("all values match expected.json", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
