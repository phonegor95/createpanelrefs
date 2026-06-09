#!/usr/bin/env python3
"""Annotate (and optionally hard-filter) deletions by loss-of-heterozygosity.

For each called deletion, counts the covered common-SNP sites (depth >= --min-dp)
and how many remain heterozygous (BAF in [het-lo, het-hi]). A real CN=1 deletion
removes one haplotype -> ~0 het sites (LOH). A depth artifact in unique sequence
keeps a normal het rate (no_LOH).

The LOH test is ONLY valid in unique, mappable sequence: inside segmental
duplications, paralogous reads create spurious het sites and fake a "no_LOH".
So any deletion overlapping the segdup track is forced to 'indeterminate_segdup'
and is NEVER hard-filtered, even with --hard-filter.

Outputs:
  --out-tsv      per-deletion LOH table
  --out-flagged  copy of the input VCF with FILTER=AllelicNoLOH set on hard-
                 filtered deletions (others untouched)
  --out-pass     hard-filtered VCF: AllelicNoLOH deletions removed (only written
                 when --hard-filter is given; otherwise identical to input pass)
"""
import argparse
import gzip


def load_bed(path):
    by = {}
    if not path:
        return by
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(("#", "@", "track", "browser")):
                continue
            p = ln.split()
            if len(p) >= 3 and p[1].isdigit():
                by.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    for c in by:
        by[c].sort()
    return by


def overlaps(idx, chrom, s, e):
    for a, b in idx.get(chrom, []):
        if a < e and b > s:
            return True
        if a >= e:
            break
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--del-bed", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--in-vcf", required=True)
    ap.add_argument("--segdup-bed", default="")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-flagged", required=True)
    ap.add_argument("--out-pass", required=True)
    ap.add_argument("--hard-filter", action="store_true")
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--min-sites", type=int, default=10)
    ap.add_argument("--het-lo", type=float, default=0.25)
    ap.add_argument("--het-hi", type=float, default=0.75)
    ap.add_argument("--loh-max-rate", type=float, default=0.02)
    ap.add_argument("--noloh-min-rate", type=float, default=0.05)
    a = ap.parse_args()

    segdup = load_bed(a.segdup_bed)

    dels = []
    with open(a.del_bed) as fh:
        for ln in fh:
            if not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            dels.append(dict(chrom=p[0], start=int(p[1]), end=int(p[2]),
                             cn=(p[3] if len(p) > 3 else "."), cov=0, het=0))

    by_chrom = {}
    with open(a.counts) as fh:
        for ln in fh:
            if ln.startswith("@") or ln.startswith("CONTIG"):
                continue
            p = ln.rstrip("\n").split("\t")
            ref, alt = int(p[2]), int(p[3])
            dp = ref + alt
            if dp < a.min_dp:
                continue
            by_chrom.setdefault(p[0], []).append((int(p[1]), alt / dp))

    # status keyed by (chrom, vcf_pos, end); vcf_pos = bed_start + 1
    status_at = {}
    for d in dels:
        for pos, baf in by_chrom.get(d["chrom"], []):
            if d["start"] < pos <= d["end"]:
                d["cov"] += 1
                if a.het_lo <= baf <= a.het_hi:
                    d["het"] += 1
        cov, het = d["cov"], d["het"]
        rate = het / cov if cov else 0.0
        if overlaps(segdup, d["chrom"], d["start"], d["end"]):
            status = "indeterminate_segdup"      # LOH test invalid in segdups
        elif cov < a.min_sites:
            status = "indeterminate_few_SNPs"
        elif rate <= a.loh_max_rate:
            status = "LOH_confirmed"
        elif rate >= a.noloh_min_rate:
            status = "no_LOH"
        else:
            status = "weak_LOH"
        d["rate"] = rate
        d["status"] = status
        status_at[(d["chrom"], d["start"] + 1, d["end"])] = status

    with open(a.out_tsv, "w") as out:
        out.write("sample\tchrom\tstart\tend\tCN\tcovered_SNPs\thet_sites\thet_rate\tLOH_status\n")
        for d in dels:
            out.write("%s\t%s\t%d\t%d\t%s\t%d\t%d\t%.3f\t%s\n" % (
                a.sample, d["chrom"], d["start"], d["end"], d["cn"],
                d["cov"], d["het"], d["rate"], d["status"]))

    FILTHDR = ('##FILTER=<ID=AllelicNoLOH,Description="CN=1 deletion in unique '
               'sequence that retains heterozygosity at common SNPs (no loss of '
               'heterozygosity) -> not allelically supported. Segdup-overlapping '
               'calls are exempt (LOH test invalid there).">\n')

    def get_end(info):
        for kv in info.split(";"):
            if kv.startswith("END="):
                return int(kv[4:])
        return None

    op = gzip.open if a.in_vcf.endswith(".gz") else open
    n_filt = 0
    with op(a.in_vcf, "rt") as fi, open(a.out_flagged, "w") as ff, open(a.out_pass, "w") as fp:
        hdr_done = False
        for line in fi:
            if line.startswith("#"):
                if line.startswith("#CHROM") and not hdr_done:
                    ff.write(FILTHDR)
                    fp.write(FILTHDR)
                    hdr_done = True
                ff.write(line)
                fp.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            is_noloh = False
            if f[4] == "<DEL>":
                end = get_end(f[7])
                if status_at.get((f[0], int(f[1]), end)) == "no_LOH":
                    is_noloh = True
            if is_noloh:
                n_filt += 1
                f[6] = "AllelicNoLOH"
                ff.write("\t".join(f) + "\n")
                if not a.hard_filter:
                    fp.write("\t".join(f) + "\n")
                # hard-filter mode: drop from the pass output
            else:
                ff.write(line)
                fp.write(line)

    print("%s: %d deletions; %d hard-filtered as no_LOH (hard_filter=%s)" % (
        a.sample, len(dels), n_filt, a.hard_filter))


if __name__ == "__main__":
    main()
