#!/usr/bin/env python3
"""Per-sample CNV QC outlier flag.

Summarises every sample's post-filter CNV profile and flags samples whose call
profile deviates from the cohort (a technical-quality red flag: noisy library,
coverage waves, GC bias). Robust statistics (median + k*MAD) are used because
the per-sample distributions are skewed.

Metrics per sample: total calls, DEL count, DUP count, CN0 (homozygous/hemizygous
loss) count, fraction of calls falling in LCR (segdup or blacklist) regions, and
total CNV bases. A sample is flagged if it exceeds median + K*MAD on total calls,
LCR calls, or CN0 calls.

This is a cohort-level QC REPORT, not a per-call filter.
Templated by Nextflow: K injected from params.
"""
import glob
import gzip
import os
from collections import defaultdict

K = float("${mad_k}")  # MAD multiplier for the outlier threshold


def load_bed(path):
    by = defaultdict(list)
    if not path or not os.path.exists(path):
        return by
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(("#", "@", "track", "browser")):
                continue
            p = ln.split()
            if len(p) >= 3 and p[1].isdigit():
                by[p[0]].append((int(p[1]), int(p[2])))
    for c in by:
        by[c].sort()
    return by


segdup = load_bed("${segdup_bed}")
blacklist = load_bed("${blacklist_bed}")


def in_region(idx, chrom, s, e):
    for a, b in idx.get(chrom, []):
        if a < e and b > s:
            return True
        if a >= e:
            break
    return False


def parse(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\\n").split("\\t")
            chrom, pos, alt, info = f[0], int(f[1]), f[4], f[7]
            end = pos
            for kv in info.split(";"):
                if kv.startswith("END="):
                    end = int(kv[4:])
            svt = "DEL" if alt == "<DEL>" else ("DUP" if alt == "<DUP>" else alt.strip("<>"))
            cn = None
            fmt = f[8].split(":")
            vals = f[9].split(":")
            if "CN" in fmt:
                cn = vals[fmt.index("CN")]
            yield chrom, pos, end, svt, cn


vcfs = sorted(glob.glob("*.filtered_segments.vcf.gz"))
per = {}
for v in vcfs:
    sid = os.path.basename(v)[:-len(".filtered_segments.vcf.gz")]
    d = dict(n=0, dele=0, dup=0, cn0=0, lcr=0, bases=0)
    for chrom, s, e, svt, cn in parse(v):
        d["n"] += 1
        d["bases"] += e - s
        if svt == "DEL":
            d["dele"] += 1
        elif svt == "DUP":
            d["dup"] += 1
        if cn == "0":
            d["cn0"] += 1
        if in_region(segdup, chrom, s, e) or in_region(blacklist, chrom, s, e):
            d["lcr"] += 1
    per[sid] = d


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return 0.0
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


def mad(xs, med):
    return median([abs(x - med) for x in xs]) or 1.0


samps = sorted(per)
thr = {}
for key in ("n", "lcr", "cn0"):
    vals = [per[s][key] for s in samps]
    med = median(vals)
    m = mad(vals, med)
    thr[key] = (med, m, med + K * m)

with open("cnv_qc_per_sample.tsv", "w") as out:
    out.write("sample\\ttotal_calls\\tDEL\\tDUP\\tCN0\\tLCR_calls\\tLCR_frac\\tCNV_Mb\\tflag\\tflag_reason\\n")
    n_flagged = 0
    for s in samps:
        d = per[s]
        reasons = []
        for key, lab in (("n", "total"), ("lcr", "LCR"), ("cn0", "CN0")):
            if d[key] > thr[key][2]:
                reasons.append("%s>%.0f" % (lab, thr[key][2]))
        flag = "OUTLIER" if reasons else "ok"
        if reasons:
            n_flagged += 1
        lcr_frac = d["lcr"] / d["n"] if d["n"] else 0.0
        out.write("%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%.2f\\t%.2f\\t%s\\t%s\\n" % (
            s, d["n"], d["dele"], d["dup"], d["cn0"], d["lcr"], lcr_frac,
            d["bases"] / 1e6, flag, ";".join(reasons) or "."))

print("CNV QC: %d samples, %d flagged as outliers (K=%.1f MAD)" % (len(samps), n_flagged, K))
for key, lab in (("n", "total_calls"), ("lcr", "LCR_calls"), ("cn0", "CN0_calls")):
    med, m, t = thr[key]
    print("  %-11s median=%.0f MAD=%.0f flag>%.0f" % (lab, med, m, t))
# Tool versions are emitted via the Nextflow `versions` topic (see main.nf).
