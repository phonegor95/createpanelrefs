#!/usr/bin/env python3
"""Cohort-contrast specificity layer for the route-A allelic CNV rescue.

The per-sample scorer (allelic_seg_call.py) is a sensitive screen: a focal run of
homozygous common SNPs in a het background flags a candidate deletion, but such
homozygous runs are common genome-wide (ROH, haplotype blocks), so per sample it
over-fires. The discriminating fact for the segdup-masked case (where depth cannot
help) is COHORT CONTRAST: a real private deletion is hemizygous in this sample at a
locus where the rest of the cohort stays HETEROZYGOUS; a systematic artifact
(segdup with no paralog-distinguishing variants, or a common ROH block) is
homozygous in ~everyone, and a common CNV recurs as a candidate in many samples.

For each per-sample del-candidate run [start, end] this computes, over the OTHER
cohort members:
  * cohort_het_support = fraction of other samples carrying >= 1 heterozygous
    common SNP inside the run. HIGH => the locus is normally heterozygous, so this
    sample's loss-of-het is anomalous (private) => a real candidate.
  * cohort_recurrence  = fraction of other samples that ALSO flag a del-candidate
    overlapping the run. HIGH => a common CNV / recurrent artifact, not private.

Classification (annotate/flag, never an auto-call):
  * private_del      : het_support >= min-het-support AND recurrence <= max-recurrence
  * recurrent_CNV    : recurrence > max-recurrence (common across the cohort)
  * uninformative_LOH: het_support < min-het-support (homozygous in ~everyone ->
                       cannot tell deletion from a no-PSV segdup / common ROH)

Inputs are the collected per-sample rescue TSVs (from allelic_seg_call.py) plus the
collected per-sample allelic-count TSVs (for the het map). Output: one cohort table
and, per sample, a contrast-annotated VCF restricted to private_del candidates.
"""
import argparse
import glob
import gzip
import os


def het_sites_from_counts(path, min_dp, het_lo, het_hi):
    """Allelic-count TSV -> {chrom: sorted [pos, ...]} of heterozygous sites."""
    by = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("@") or ln.startswith("CONTIG"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            ref, alt = int(f[2]), int(f[3])
            dp = ref + alt
            if dp < min_dp:
                continue
            if het_lo <= alt / dp <= het_hi:
                by.setdefault(f[0], []).append(int(f[1]))
    for c in by:
        by[c].sort()
    return by


def has_het_in(hets, chrom, s, e):
    """True if sample `hets` has any heterozygous site in (s, e]."""
    import bisect
    lst = hets.get(chrom)
    if not lst:
        return False
    i = bisect.bisect_right(lst, s)
    return i < len(lst) and lst[i] <= e


def load_candidates(path, classes):
    """rescue TSV -> list of dict del-candidate rows (header-driven)."""
    rows = []
    with open(path) as fh:
        hdr = None
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if hdr is None:
                hdr = {n: i for i, n in enumerate(f)}
                continue
            if f[hdr["class"]] not in classes:
                continue
            rows.append(dict(chrom=f[hdr["chrom"]], start=int(f[hdr["start"]]),
                             end=int(f[hdr["end"]]), cls=f[hdr["class"]],
                             conf=f[hdr["confidence"]], line=f, hdr=hdr))
    return rows


def overlap(a_s, a_e, b_s, b_e):
    return a_s < b_e and b_s < a_e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rescue-glob", required=True,
                    help="glob for per-sample *.rescue.tsv (e.g. 'rescue_out/*.rescue.tsv')")
    ap.add_argument("--counts-glob", required=True,
                    help="glob for per-sample allelic-count TSVs; basename must start "
                         "with the sample id (e.g. 'ac.<sample>.tsv')")
    ap.add_argument("--counts-prefix", default="ac.")
    ap.add_argument("--counts-suffix", default=".tsv")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-dir", default="rescue_contrast")
    ap.add_argument("--min-dp", type=int, default=8)
    ap.add_argument("--het-lo", type=float, default=0.25)
    ap.add_argument("--het-hi", type=float, default=0.75)
    ap.add_argument("--min-het-support", type=float, default=0.30,
                    help="fraction of OTHER samples that must stay het over the run")
    ap.add_argument("--max-recurrence", type=float, default=0.80,
                    help="fraction of OTHER samples flagging the run above which it is "
                         "a common CNV / recurrent artifact. Kept high so a carrier-"
                         "enriched locus in a SMALL cohort is not wrongly suppressed; "
                         "in a large screening cohort recurrence for a rare variant is "
                         "low and this only catches near-universal common CNVs.")
    ap.add_argument("--del-classes", default="depth_silent_del_candidate,"
                    "segdup_masked_del_candidate,gcnv_concordant_del")
    a = ap.parse_args()

    classes = set(a.del_classes.split(","))

    # sample id <- rescue file (first column of the body is the sample)
    rescue_files = sorted(glob.glob(a.rescue_glob))
    samples = {}
    for rf in rescue_files:
        with open(rf) as fh:
            next(fh, None)
            first = next(fh, None)
        sid = first.split("\t")[0] if first else os.path.basename(rf).split(".")[0]
        samples[sid] = dict(rescue=rf)

    # match each sample to its allelic-count file by id substring
    for cf in glob.glob(a.counts_glob):
        base = os.path.basename(cf)
        sid = base[len(a.counts_prefix):-len(a.counts_suffix)] \
            if base.startswith(a.counts_prefix) and base.endswith(a.counts_suffix) else None
        if sid in samples:
            samples[sid]["counts"] = cf

    het_maps = {sid: het_sites_from_counts(s["counts"], a.min_dp, a.het_lo, a.het_hi)
                for sid, s in samples.items() if "counts" in s}
    cand = {sid: load_candidates(s["rescue"], classes) for sid, s in samples.items()}

    n_other = max(len(samples) - 1, 1)
    os.makedirs(a.out_dir, exist_ok=True)
    out = open(a.out_tsv, "w")
    out.write("sample\tchrom\tstart\tend\tclass\tconf\tcohort_het_support"
              "\tcohort_recurrence\tcontrast_class\n")

    summary = {}
    for sid in samples:
        priv = open(os.path.join(a.out_dir, sid + ".rescue.private.vcf"), "w")
        priv.write("##fileformat=VCFv4.2\n##ALT=<ID=DEL,Description=\"Deletion\">\n")
        priv.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n')
        priv.write('##INFO=<ID=COHORTHETSUPPORT,Number=1,Type=Float,Description="fraction of other cohort samples het over the run">\n')
        priv.write('##INFO=<ID=COHORTRECURRENCE,Number=1,Type=Float,Description="fraction of other cohort samples flagging the run">\n')
        priv.write('##INFO=<ID=CONTRASTCLASS,Number=1,Type=String,Description="cohort-contrast class">\n')
        priv.write('##FILTER=<ID=AllelicRescuePrivate,Description="Private hemizygous focal LOH: deletion candidate the depth caller missed, anomalous vs the cohort.">\n')
        priv.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        cnt = dict(private_del=0, recurrent_CNV=0, uninformative_LOH=0)
        for c in cand[sid]:
            het_sup = sum(1 for o in samples if o != sid
                          and has_het_in(het_maps.get(o, {}), c["chrom"], c["start"], c["end"]))
            recur = sum(1 for o in samples if o != sid
                        and any(overlap(c["start"], c["end"], x["start"], x["end"])
                                for x in cand.get(o, [])))
            hs, rc = het_sup / n_other, recur / n_other
            if rc > a.max_recurrence:
                cc = "recurrent_CNV"
            elif hs >= a.min_het_support:
                cc = "private_del"
            else:
                cc = "uninformative_LOH"
            cnt[cc] += 1
            out.write("%s\t%s\t%d\t%d\t%s\t%s\t%.3f\t%.3f\t%s\n" % (
                sid, c["chrom"], c["start"], c["end"], c["cls"], c["conf"], hs, rc, cc))
            if cc == "private_del":
                priv.write("%s\t%d\t.\tN\t<DEL>\t.\tAllelicRescuePrivate\t"
                           "END=%d;COHORTHETSUPPORT=%.3f;COHORTRECURRENCE=%.3f;"
                           "CONTRASTCLASS=%s\n"
                           % (c["chrom"], c["start"], c["end"], hs, rc, cc))
        priv.close()
        summary[sid] = cnt
    out.close()

    print("%-22s %8s %10s %14s" % ("sample", "private", "recurrent", "uninformative"))
    for sid, cnt in summary.items():
        print("%-22s %8d %10d %14d"
              % (sid, cnt["private_del"], cnt["recurrent_CNV"], cnt["uninformative_LOH"]))


if __name__ == "__main__":
    main()
