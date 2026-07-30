#!/usr/bin/env python3
"""Cohort-recurrence (two-arm) soft filter for gCNV genotyped segments.

Flags calls that recur across the cohort as likely artifacts / common CNVs,
without dropping any record (soft FILTER=CohortRecurrent). Two arms:

  * reciprocal-coverage arm: a call is "shared" by another sample when the
    overlap covers >= COV of THIS call; flag if shared-fraction >= T1.
  * any-overlap arm: a call is "touched" by another sample on any base; flag
    if touched-fraction >= T2. Catches ragged-breakpoint paralog artifacts
    (e.g. GTF2I) that escape the reciprocal arm.

Grouping: autosomes are pooled over the whole cohort (sex-agnostic); chrX and
chrY are split by inferred sex, because hemizygous male allosome calls would
otherwise look "recurrent" against the female cohort and vice-versa.

Denominators are derived from the cohort sex map at run time (no hard-coding).
Templated by Nextflow: T1/COV/T2 are injected from params.
"""
import glob
import gzip
import os
import subprocess
from collections import defaultdict

# --- thresholds injected by Nextflow template engine -----------------------
T1  = float("${recip_threshold}")    # reciprocal-coverage arm
COV = float("${recip_coverage}")     # coverage fraction defining "shared"
T2  = float("${overlap_threshold}")  # any-overlap arm
MIN_PEERS = int("${min_cohort}")     # min peers per call for a meaningful frequency
SEXCHR = {"chrX", "chrY", "X", "Y"}

OUTDIR = "recurfilt"
os.makedirs(OUTDIR, exist_ok=True)

# --- cohort sex map ---------------------------------------------------------
sex_of = {}
for ln in open("${sex_map}"):
    ln = ln.rstrip("\\n")
    if not ln:
        continue
    sid, sex = ln.split("\\t")[:2]
    sex_of[sid] = sex
NALL = len(sex_of)
NSEX = {"F": sum(1 for s in sex_of.values() if s == "F"),
        "M": sum(1 for s in sex_of.values() if s == "M")}
print("cohort: %d samples (%dF / %dM)" % (NALL, NSEX["F"], NSEX["M"]))

FILTLINE = (
    '##FILTER=<ID=CohortRecurrent,Description="Cohort-recurrent CNV/artifact '
    '(two-arm): reciprocal-coverage(>=%.2f) freq >= %.2f, OR any-overlap freq '
    '>= %.2f. Frequencies are over the PEER samples (the carrier itself is '
    'excluded from numerator and denominator). Autosomes pooled over the whole '
    'cohort; chrX/chrY split by inferred sex.">\\n' % (COV, T1, T2)
)


def parse_vcf(path):
    """Yield (chrom, start, end, svt, raw_line_fields) for each call."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fi:
        for line in fi:
            if line.startswith("#"):
                continue
            f = line.rstrip("\\n").split("\\t")
            chrom, pos, alt, info = f[0], int(f[1]), f[4], f[7]
            if alt == "<DEL>":
                svt = "DEL"
            elif alt == "<DUP>":
                svt = "DUP"
            else:
                svt = alt.strip("<>")
            end = None
            for kv in info.split(";"):
                if kv.startswith("END="):
                    end = int(kv[4:])
            if end is None:
                end = pos
            yield chrom, pos, end, svt, f


# --- gather every sample's calls -------------------------------------------
vcfs = sorted(glob.glob("*.filtered_segments.vcf.gz"))
calls = []  # dict per call
for v in vcfs:
    sid = os.path.basename(v)[:-len(".filtered_segments.vcf.gz")]
    if sid not in sex_of:
        raise SystemExit("sample %r in VCFs but not in sex map" % sid)
    for chrom, start, end, svt, _ in parse_vcf(v):
        calls.append(dict(chrom=chrom, start=start, end=end, svt=svt,
                          samp=sid, sex=sex_of[sid]))
print("loaded %d calls from %d samples" % (len(calls), len(vcfs)))


def grpkey(c):
    if c["chrom"] in SEXCHR:
        return (c["sex"], c["chrom"], c["svt"]), NSEX[c["sex"]]
    return ("ALL", c["chrom"], c["svt"]), NALL


grp = defaultdict(list)
for c in calls:
    grp[grpkey(c)[0]].append(c)
for k in grp:
    grp[k].sort(key=lambda c: c["start"])


def freqs(c):
    """Fraction of the OTHER samples in this call's group that share / touch it.

    The call's own sample is excluded from both numerator and denominator. With
    self-counting a private call always scores 1/denom, so at the default
    T1=0.5 a 2-sample cohort would flag every single call, and a lone female in
    a male cohort would have every chrX call flagged against a denominator of 1.
    """
    key, denom = grpkey(c)
    peers = denom - 1
    if peers < 1:
        return 0.0, 0.0
    g = grp[key]
    need = COV * (c["end"] - c["start"])
    rec, ovl = set(), set()
    for o in g:
        if o["start"] >= c["end"]:
            break
        if o["end"] <= c["start"] or o["samp"] == c["samp"]:
            continue
        inter = min(o["end"], c["end"]) - max(o["start"], c["start"])
        ovl.add(o["samp"])
        if inter >= need:
            rec.add(o["samp"])
    return len(rec) / peers, len(ovl) / peers


# A recurrence frequency over a handful of peers is quantised too coarsely to
# mean anything (with 3 peers the only reachable values are 0, .33, .67, 1).
# Rather than flag near-arbitrarily, degrade to a no-op and say so loudly.
flag = set()
if NALL - 1 < MIN_PEERS:
    print("WARNING: %d-sample cohort leaves only %d peer(s) per call, below the "
          "minimum of %d (params.recur_min_cohort). The recurrence frequency is "
          "not meaningful at this size: NO calls will be flagged CohortRecurrent "
          "and every call is emitted as PASS." % (NALL, NALL - 1, MIN_PEERS))
else:
    for c in calls:
        rc, ao = freqs(c)
        if rc >= T1 or ao >= T2:
            flag.add((c["samp"], c["chrom"], c["start"]))
    print("flagged %d / %d calls as CohortRecurrent" % (len(flag), len(calls)))

# --- write per-sample soft-flagged + pass-only VCFs ------------------------
stats = []
for v in vcfs:
    sid = os.path.basename(v)[:-len(".filtered_segments.vcf.gz")]
    soft = os.path.join(OUTDIR, sid + ".recurfilt.vcf")
    passf = os.path.join(OUTDIR, sid + ".recurfilt.pass.vcf")
    nflag = ntot = npass = 0
    hdr_done = False
    with gzip.open(v, "rt") as fi, open(soft, "w") as fo, open(passf, "w") as fp:
        for line in fi:
            if line.startswith("#"):
                if line.startswith("#CHROM") and not hdr_done:
                    fo.write(FILTLINE)
                    fp.write(FILTLINE)
                    hdr_done = True
                fo.write(line)
                fp.write(line)
                continue
            f = line.rstrip("\\n").split("\\t")
            ntot += 1
            pos = int(f[1])
            if (sid, f[0], pos) in flag:
                f[6] = "CohortRecurrent"
                nflag += 1
            else:
                f[6] = "PASS" if f[6] in (".", "") else f[6]
            fo.write("\\t".join(f) + "\\n")
            # PASS-only output: emit genuinely-passing records, not merely
            # "not CohortRecurrent" (a pre-existing non-PASS FILTER must not leak).
            if f[6] == "PASS":
                fp.write("\\t".join(f) + "\\n")
                npass += 1
    for path in (soft, passf):
        subprocess.run(["bgzip", "-f", path], check=True)
        subprocess.run(["tabix", "-p", "vcf", "-f", path + ".gz"], check=True)
    stats.append((sid, ntot, nflag, npass))

print("%-34s %6s %8s %6s" % ("sample", "total", "flagged", "PASS"))
tt = tf = tp = 0
for sid, ntot, nflag, npass in stats:
    print("%-34s %6d %8d %6d" % (sid, ntot, nflag, npass))
    tt += ntot
    tf += nflag
    tp += npass
print("%-34s %6d %8d %6d" % ("TOTAL", tt, tf, tp))
# Tool versions are emitted via the Nextflow `versions` topic (see main.nf).
