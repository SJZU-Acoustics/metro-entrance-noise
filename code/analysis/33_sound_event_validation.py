#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
33_sound_event_validation.py

Independently recompute the manual-validation metrics for the automated
sound-event classifier (Paper 1 Methods), from the per-clip
validation table data/validation/sound_event_validation_clips.csv

Reference standard = majority vote of three blind annotators (A, B, C);
blanks treated as 'unknown'; a clip with no >=2 agreement is excluded from
that metric (the first author's protocol, docx D.3).

Two resolutions are computed:
  (1) 12-category event code  -> verifies the first author's reported numbers
  (2) 4 source groups (traffic / human / natural / mechanical)
      -> the directly paper-relevant reliability (station-level source shares
         are built from the 4 groups, so within-group confusions are not errors)

For dominant and secondary source, overall and for the night subsample, we
report N (valid-majority clips), accuracy, Cohen's kappa, the model 'unknown'
share, and per-group precision/recall/F1 (4-group, dominant).

Outputs (exploration/output/):
  33_sound_event_validation__summary.csv      headline metrics
  33_sound_event_validation__perclass.csv     4-group per-class P/R/F1 (dominant)
  33_sound_event_validation__confusion.csv    4-group confusion (dominant, overall)
Plus a human-readable brief: exploration/scripts/../reports? -> printed to stdout
"""
import os
from collections import Counter
from sklearn.metrics import (accuracy_score, cohen_kappa_score,
                             precision_recall_fscore_support, confusion_matrix)
import csv

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repository root
CLIPS_CSV = os.path.join(ROOT, "data", "validation", "sound_event_validation_clips.csv")
OUTDIR = os.path.join(ROOT, "intermediate", "analysis")
os.makedirs(OUTDIR, exist_ok=True)

# 12-category code -> 4 source group (from docx D.1 / main.tex Methods)
CODE2GROUP = {
    1: "natural", 2: "natural",                    # birdsong, insect
    3: "traffic", 4: "traffic", 5: "traffic",      # car pass-by, horn, motorcycle
    6: "human", 7: "human", 8: "human", 9: "human", 11: "human",
    #          conversation, music, vendor call, shop loudspeaker, footsteps
    10: "mechanical", 12: "mechanical",            # cleaning, mechanical friction
}
UNK = "unknown"


def norm_code(v):
    """Cell -> int code, or UNK for blank."""
    if v is None:
        return UNK
    if isinstance(v, str):
        v = v.strip()
        if v == "":
            return UNK
        try:
            return int(float(v))
        except ValueError:
            return v
    if isinstance(v, float):
        return int(v)
    return v


def to_group(code):
    if code == UNK:
        return UNK
    return CODE2GROUP.get(code, UNK)


def majority(labels):
    """Majority among 3 labels; None if no value reaches count>=2."""
    c = Counter(labels)
    val, n = c.most_common(1)[0]
    return val if n >= 2 else None


# ---- load raw clips ----
with open(CLIPS_CSV, newline="", encoding="utf-8") as f:
    rdr = csv.reader(f)
    hdr = next(rdr)
    rows = list(rdr)
idx = {name: i for i, name in enumerate(hdr)}

clips = []
for r in rows:
    if not r[idx["clip_id"]]:
        continue
    clips.append({
        "clip_id": r[idx["clip_id"]],
        "station": str(r[idx["source_audio_file"]]).split("_")[0] + "_" +
                   str(r[idx["source_audio_file"]]).split("_")[1],
        "day_night": (r[idx["day_night"]] or "").strip(),
        "m_pri": norm_code(r[idx["model_primary_event_code"]]),
        "a_pri": norm_code(r[idx["A_primary_event"]]),
        "b_pri": norm_code(r[idx["B_primary_event"]]),
        "c_pri": norm_code(r[idx["C_primary_event"]]),
        "m_sec": norm_code(r[idx["model_secondary_event_code"]]),
        "a_sec": norm_code(r[idx["A_secondary_event"]]),
        "b_sec": norm_code(r[idx["B_secondary_event"]]),
        "c_sec": norm_code(r[idx["C_secondary_event"]]),
    })

print(f"Loaded {len(clips)} clips")
stations = Counter(c["station"] for c in clips)
print("Source stations:", dict(stations))
dn = Counter(c["day_night"] for c in clips)
print("day/night split:", dict(dn))


def evaluate(clips, task, level):
    """
    task: 'pri' or 'sec'; level: 'code' (12-cat) or 'group' (4-group).
    Returns dict of metrics over clips with a valid annotator majority.
    """
    mk = "m_" + task
    ak, bk, ck = "a_" + task, "b_" + task, "c_" + task
    y_true, y_pred = [], []
    model_unknown = 0
    n_total = len(clips)
    for c in clips:
        if level == "group":
            anns = [to_group(c[ak]), to_group(c[bk]), to_group(c[ck])]
            mdl = to_group(c[mk])
        else:
            anns = [c[ak], c[bk], c[ck]]
            mdl = c[mk]
        maj = majority(anns)
        if maj is None:
            continue
        # stringify so sklearn never sees a str/int label mix (codes + 'unknown')
        y_true.append(str(maj))
        y_pred.append(str(mdl))
        if mdl == UNK:
            model_unknown += 1
    n = len(y_true)
    acc = accuracy_score(y_true, y_pred) if n else float("nan")
    kappa = cohen_kappa_score(y_true, y_pred) if n else float("nan")
    return {
        "task": task, "level": level,
        "N_total": n_total, "N_valid_majority": n,
        "accuracy": round(acc, 3), "kappa": round(kappa, 3),
        "model_unknown_share_of_valid": round(model_unknown / n, 3) if n else float("nan"),
        "model_unknown_share_of_total": round(
            sum(1 for c in clips if (to_group(c[mk]) if level == "group" else c[mk]) == UNK)
            / n_total, 3),
        "_y_true": y_true, "_y_pred": y_pred,
    }


def show(d):
    print(f"  [{d['level']:5s} | {d['task']}] N_valid={d['N_valid_majority']}/{d['N_total']}  "
          f"acc={d['accuracy']}  kappa={d['kappa']}  "
          f"model_unk(of valid)={d['model_unknown_share_of_valid']}  "
          f"model_unk(of total)={d['model_unknown_share_of_total']}")


night = [c for c in clips if c["day_night"] == "night"]
day = [c for c in clips if c["day_night"] == "day"]

summary_rows = []
print("\n=== OVERALL (all 300 clips) ===")
for level in ("code", "group"):
    for task in ("pri", "sec"):
        d = evaluate(clips, task, level)
        show(d)
        summary_rows.append({k: v for k, v in d.items() if not k.startswith("_")} | {"subset": "all"})

print("\n=== NIGHT subsample ===")
for level in ("code", "group"):
    for task in ("pri", "sec"):
        d = evaluate(night, task, level)
        show(d)
        summary_rows.append({k: v for k, v in d.items() if not k.startswith("_")} | {"subset": "night"})

print("\n=== DAY subsample (for completeness) ===")
for level in ("code", "group"):
    for task in ("pri", "sec"):
        d = evaluate(day, task, level)
        show(d)
        summary_rows.append({k: v for k, v in d.items() if not k.startswith("_")} | {"subset": "day"})

# ---- per-class P/R/F1, 4-group dominant, overall ----
d = evaluate(clips, "pri", "group")
labels = ["traffic", "human", "natural", "mechanical", UNK]
labels = [l for l in labels if l in set(d["_y_true"]) | set(d["_y_pred"])]
P, R, F, S = precision_recall_fscore_support(d["_y_true"], d["_y_pred"],
                                             labels=labels, zero_division=0)
perclass = []
for l, p, r, f, s in zip(labels, P, R, F, S):
    perclass.append({"group": l, "precision": round(p, 3), "recall": round(r, 3),
                     "f1": round(f, 3), "support_majority": int(s)})
print("\n=== Per-class P/R/F1 (4-group, dominant, overall) ===")
for pc in perclass:
    print(f"  {pc['group']:11s}  P={pc['precision']:.3f}  R={pc['recall']:.3f}  "
          f"F1={pc['f1']:.3f}  support={pc['support_majority']}")

cm = confusion_matrix(d["_y_true"], d["_y_pred"], labels=labels)
print("\n=== Confusion (rows=human majority, cols=model), 4-group dominant overall ===")
print("            " + "  ".join(f"{l[:5]:>5s}" for l in labels))
for l, row in zip(labels, cm):
    print(f"  {l:11s}" + "  ".join(f"{v:5d}" for v in row))

# ---- inter-annotator agreement (sets the human ceiling the model is judged against) ----
def pairwise_kappa(clips, task, level, k1, k2):
    y1, y2 = [], []
    for c in clips:
        v1 = to_group(c[k1 + "_" + task]) if level == "group" else c[k1 + "_" + task]
        v2 = to_group(c[k2 + "_" + task]) if level == "group" else c[k2 + "_" + task]
        y1.append(str(v1)); y2.append(str(v2))
    return round(cohen_kappa_score(y1, y2), 3)

print("\n=== Inter-annotator Cohen's kappa (dominant source) ===")
inter_rows = []
for level in ("code", "group"):
    ks = []
    for a, b in (("a", "b"), ("a", "c"), ("b", "c")):
        k = pairwise_kappa(clips, "pri", level, a, b)
        ks.append(k)
        inter_rows.append({"level": level, "pair": f"{a.upper()}-{b.upper()}", "kappa": k})
    mean_k = round(sum(ks) / len(ks), 3)
    inter_rows.append({"level": level, "pair": "mean_human", "kappa": mean_k})
    # model vs each human
    mks = [pairwise_kappa(clips, "pri", level, "m", h) for h in ("a", "b", "c")]
    for h, mk in zip(("a", "b", "c"), mks):
        inter_rows.append({"level": level, "pair": f"M-{h.upper()}", "kappa": mk})
    print(f"  [{level:5s}] human pairwise A-B/A-C/B-C = {ks}, mean = {mean_k}; "
          f"model-vs-human M-A/M-B/M-C = {mks}, mean = {round(sum(mks)/3,3)}")

inter_path = os.path.join(OUTDIR, "33_sound_event_validation__interannotator.csv")
with open(inter_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["level", "pair", "kappa"])
    w.writeheader()
    w.writerows(inter_rows)

# ---- AGGREGATE-SHARE / BIAS analysis (the paper-relevant unit) ----
# The manuscript uses station-level 4-group *shares* built from the dominant
# event, not per-clip labels. What matters for the headline (82% traffic) is
# whether the model's aggregate shares are *unbiased* vs the human consensus.
# Restrict to clips with a human-majority in one of the 4 groups (drop the few
# unknown-majority clips) so model and majority share-vectors both sum to 1.
import numpy as np

GROUPS = ["traffic", "human", "natural", "mechanical"]


def share_vec(labels):
    n = len(labels)
    return {g: sum(1 for x in labels if x == g) / n for g in GROUPS}


paired = []  # (model_group, majority_group) on clips with majority in 4 groups
for c in clips:
    anns = [to_group(c["a_pri"]), to_group(c["b_pri"]), to_group(c["c_pri"])]
    maj = majority(anns)
    if maj in GROUPS:
        paired.append((to_group(c["m_pri"]), maj))

m_lab = [p[0] for p in paired]
h_lab = [p[1] for p in paired]
n_pair = len(paired)
m_share = share_vec(m_lab)
h_share = share_vec(h_lab)

# per-annotator shares (on the same clip set) to show the model sits in human spread
ann_shares = {}
for who, key in (("A", "a_pri"), ("B", "b_pri"), ("C", "c_pri")):
    labs = []
    for c in clips:
        anns = [to_group(c["a_pri"]), to_group(c["b_pri"]), to_group(c["c_pri"])]
        if majority(anns) in GROUPS:
            labs.append(to_group(c[key]))
    ann_shares[who] = share_vec(labs)

# bootstrap 1000x (house rule) -> 95% CI on (model - majority) share difference
rng = np.random.default_rng(11)
B = 1000
idx = np.arange(n_pair)
boot_diff = {g: [] for g in GROUPS}
for _ in range(B):
    samp = rng.choice(idx, size=n_pair, replace=True)
    ms = share_vec([m_lab[i] for i in samp])
    hs = share_vec([h_lab[i] for i in samp])
    for g in GROUPS:
        boot_diff[g].append(ms[g] - hs[g])

print(f"\n=== AGGREGATE 4-group SHARE agreement (dominant; {n_pair} clips with a 4-group human majority) ===")
print(f"  {'group':11s} {'model%':>7s} {'human%':>7s} {'diff(pp)':>9s} {'95% CI (pp)':>18s}   annotators A/B/C %")
share_rows = []
for g in GROUPS:
    d = (m_share[g] - h_share[g]) * 100
    lo, hi = np.percentile(boot_diff[g], [2.5, 97.5]) * 100
    abc = "/".join(f"{ann_shares[w][g]*100:.0f}" for w in ("A", "B", "C"))
    print(f"  {g:11s} {m_share[g]*100:7.1f} {h_share[g]*100:7.1f} {d:+9.1f} "
          f"   [{lo:+.1f}, {hi:+.1f}]      {abc}")
    share_rows.append({"group": g, "model_share_pct": round(m_share[g]*100, 1),
                       "human_majority_share_pct": round(h_share[g]*100, 1),
                       "diff_pp": round(d, 2), "ci_lo_pp": round(lo, 2),
                       "ci_hi_pp": round(hi, 2),
                       "ci_includes_zero": bool(lo <= 0 <= hi),
                       "annA_pct": round(ann_shares['A'][g]*100, 1),
                       "annB_pct": round(ann_shares['B'][g]*100, 1),
                       "annC_pct": round(ann_shares['C'][g]*100, 1)})
tvd = 0.5 * sum(abs(m_share[g] - h_share[g]) for g in GROUPS)
print(f"  Total variation distance (model vs human composition) = {tvd:.3f}  "
      f"(0 = identical; max 1)")
print("  NB: validation sample is balanced across groups, NOT citywide-representative "
      f"(traffic {h_share['traffic']*100:.0f}% here vs 82% citywide).")

sh_path = os.path.join(OUTDIR, "33_sound_event_validation__shareagreement.csv")
with open(sh_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(share_rows[0].keys()))
    w.writeheader(); w.writerows(share_rows)

# ---- write outputs ----
sum_path = os.path.join(OUTDIR, "33_sound_event_validation__summary.csv")
with open(sum_path, "w", newline="") as f:
    cols = ["subset", "task", "level", "N_total", "N_valid_majority",
            "accuracy", "kappa", "model_unknown_share_of_valid",
            "model_unknown_share_of_total"]
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for r in summary_rows:
        w.writerow({k: r[k] for k in cols})

pc_path = os.path.join(OUTDIR, "33_sound_event_validation__perclass.csv")
with open(pc_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["group", "precision", "recall", "f1", "support_majority"])
    w.writeheader()
    w.writerows(perclass)

cm_path = os.path.join(OUTDIR, "33_sound_event_validation__confusion.csv")
with open(cm_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["human_majority\\model"] + labels)
    for l, row in zip(labels, cm):
        w.writerow([l] + list(row))

print(f"\nWrote:\n  {sum_path}\n  {pc_path}\n  {cm_path}")
