#!/usr/bin/env python3
"""Annotate per-sample deletions with loss-of-heterozygosity (LOH) evidence.

Reads a BED of called deletions and a GATK CollectAllelicCounts table, then for
each deletion counts the common-SNP sites it covers (depth >= --min-dp) and how
many are heterozygous (B-allele fraction in [0.25, 0.75]). A real CN=1 deletion
removes one haplotype -> the remaining haplotype is the only signal -> ~0 het
sites (LOH). This is an orthogonal, allele-based confirmation independent of the
read-depth signal gCNV used. Only meaningful in unique, mappable sequence; in
segdups/LCRs het genotypes are unreliable, so we mark those 'indeterminate'.

Soft annotation only: never drops a call.
"""
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--del-bed", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--min-sites", type=int, default=10,
                    help="min covered SNPs for an LOH/no-LOH verdict (else indeterminate)")
    ap.add_argument("--het-lo", type=float, default=0.25)
    ap.add_argument("--het-hi", type=float, default=0.75)
    ap.add_argument("--loh-max-rate", type=float, default=0.02,
                    help="het-rate <= this => LOH_confirmed")
    ap.add_argument("--noloh-min-rate", type=float, default=0.05,
                    help="het-rate >= this => no_LOH")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    dels = []
    with open(a.del_bed) as fh:
        for ln in fh:
            if not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            dels.append(dict(chrom=p[0], start=int(p[1]), end=int(p[2]),
                             cn=(p[3] if len(p) > 3 else "."), cov=0, het=0))

    # bucket counts by chrom for speed
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

    for d in dels:
        for pos, baf in by_chrom.get(d["chrom"], []):
            if d["start"] < pos <= d["end"]:
                d["cov"] += 1
                if a.het_lo <= baf <= a.het_hi:
                    d["het"] += 1

    with open(a.out, "w") as out:
        out.write("sample\tchrom\tstart\tend\tCN\tcovered_SNPs\thet_sites\thet_rate\tLOH_status\n")
        for d in dels:
            cov, het = d["cov"], d["het"]
            rate = het / cov if cov else 0.0
            if cov < a.min_sites:
                status = "indeterminate_few_SNPs"
            elif rate <= a.loh_max_rate:
                status = "LOH_confirmed"
            elif rate >= a.noloh_min_rate:
                status = "no_LOH"
            else:
                status = "weak_LOH"
            out.write("%s\t%s\t%d\t%d\t%s\t%d\t%d\t%.3f\t%s\n" % (
                a.sample, d["chrom"], d["start"], d["end"], d["cn"], cov, het, rate, status))


if __name__ == "__main__":
    main()
