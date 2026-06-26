#!/usr/bin/env python3
"""Route-A allelic/zygosity rescue scanner for segdup-masked focal deletions.

Depth-based gCNV (and ModelSegments) are structurally blind to deletions
embedded in segmental duplications -- the canonical case being the alpha-globin
-a3.7 / -a4.2 thalassemia deletions. The 1 kb bin averages in paralog-mapped
reads so the half-copy dip never appears, and ModelSegments both denoises the CR
axis flat AND drops the hemizygous SNPs from its MAF fit, so the locus is swept
into the flanking diploid arm. The per-SNP allelic axis is the only instrument
that keeps the signal, and that is what this scanner reads directly.

Two facts hold at uniquely-mappable common-SNP positions when CollectAllelicCounts
is run at high min-MAPQ (>= 30) on the recalibrated CRAM:

  * a heterozygous-background deletion forces HEMIZYGOSITY -> the site reads as
    homozygous (BAF -> 0 or 1) even though the germline genotype is het; and
  * it HALVES the unique-read depth at that site.

So a *focal run* of common SNPs that are simultaneously homozygous-looking AND
depth-reduced, sitting in an otherwise heterozygous background, is a deletion.

Two refinements over a naive global-median scan (which over-fired badly: a single
genome-wide depth threshold flags every GC/mappability coverage trough):

  1. LOCAL ROLLING-MEDIAN DEPTH BASELINE. Each site's depth is judged against the
     median depth of its local neighbourhood (a ~100 kb tile), NOT the genome-wide
     median. Regional coverage waves are absorbed into the baseline, so a deletion
     must be a *local* dip relative to its own neighbours -- the only thing a true
     hemizygous loss actually is.

  2. LOH-GATING. Loss of heterozygosity is necessary but NOT sufficient: copy-
     neutral LOH and runs-of-homozygosity (ROH) are hemizygous-looking too, but at
     NORMAL depth. The depth gate already separates those (deletion ratio ~0.5,
     ROH ~1.0). LOH-gating adds the converse guard: the run must sit in a genuinely
     HETEROZYGOUS background (its flanks must carry het sites at a normal rate),
     otherwise the loss-of-het inside the run is uninformative (a long ROH block,
     not a focal deletion). Deletion = depth-reduced + hemizygous + het-flanked.

Output is an annotation/flag for review (a candidate-DEL TSV and an optional VCF),
never a silent auto-call -- per design, the depth-based caller stays authoritative
and this rescues what it cannot see.
"""
import argparse
import gzip
import os


# --------------------------------------------------------------------------- IO
def parse_allelic_counts(path, min_dp):
    """GATK CollectAllelicCounts TSV -> {chrom: [(pos, ref, alt), ...]} sorted."""
    by = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("@") or ln.startswith("CONTIG"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            chrom, pos, ref, alt = f[0], int(f[1]), int(f[2]), int(f[3])
            if ref + alt < min_dp:
                continue
            by.setdefault(chrom, []).append((pos, ref, alt))
    for c in by:
        by[c].sort()
    return by


def load_bed(path):
    """BED -> {chrom: [(start, end), ...]} sorted (for the segdup track)."""
    by = {}
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
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


def overlaps_bed(idx, chrom, s, e):
    for a, b in idx.get(chrom, []):
        if a < e and b > s:
            return True
        if a >= e:
            break
    return False


def parse_gcnv_segments(path):
    """gCNV genotyped_segments VCF -> {chrom: [(start, end, cn, svt), ...]}."""
    by = {}
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return by
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            chrom, pos, alt, info = f[0], int(f[1]), f[4], f[7]
            end = pos
            for kv in info.split(";"):
                if kv.startswith("END="):
                    end = int(kv[4:])
            cn = None
            fmt, val = f[8].split(":"), f[9].split(":")
            if "CN" in fmt:
                try:
                    cn = int(val[fmt.index("CN")])
                except ValueError:
                    cn = None
            svt = "DEL" if alt == "<DEL>" else ("DUP" if alt == "<DUP>" else ".")
            by.setdefault(chrom, []).append((pos, end, cn, svt))
    for c in by:
        by[c].sort()
    return by


def gcnv_cn_over(gcnv, chrom, s, e):
    """Minimum gCNV CN of any segment overlapping [s, e] (None if no coverage)."""
    best = None
    for ps, pe, cn, _svt in gcnv.get(chrom, []):
        if ps > e:
            break
        if pe < s or cn is None:
            continue
        best = cn if best is None else min(best, cn)
    return best


# ----------------------------------------------------------------- math helpers
def median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return 0.0
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


def local_baselines(positions, depths, window_bp):
    """Local rolling-median depth baseline per site (refinement #1).

    Tile the chromosome into non-overlapping `window_bp` windows and take the
    median depth of each tile; a focal deletion perturbs only a handful of the
    ~hundreds of sites in its tile, so the tile median tracks regional coverage,
    not the deletion. Each site inherits the median of the tile it falls in, with
    a 3-tile (prev/this/next) pool so the baseline is smooth across boundaries and
    robust where a tile is sparse.
    """
    tiles = {}
    for pos, dp in zip(positions, depths):
        tiles.setdefault(pos // window_bp, []).append(dp)
    tile_med = {}
    for t in tiles:
        pool = tiles.get(t - 1, []) + tiles[t] + tiles.get(t + 1, [])
        tile_med[t] = median(pool) or 1.0
    return [tile_med[pos // window_bp] for pos in positions]


# ----------------------------------------------------------------- core scanning
def scan_chrom(chrom, lst, a, gcnv, segdup):
    """Yield candidate loss-of-heterozygosity rows for one chromosome.

    PRIMARY signal = a run of consecutive HOMOZYGOUS common SNPs (loss of het),
    detected depth-independently. This is the robust axis: a hemizygous deletion
    forces every covered SNP to read homozygous.

    Depth behaves OPPOSITELY in unique sequence vs segdups, so it cannot be the
    primary discriminant -- only a classifier:
      * in UNIQUE sequence a hemizygous deletion HALVES depth; and
      * in a SEGMENTAL DUPLICATION the ~95%-identical surviving paralog keeps
        mapping reads to the deleted gene's coordinates, so depth stays ~NORMAL
        even though zygosity collapses to hemizygous (the alpha-globin case).

    GATING (the two refinements):
      * local rolling-median DEPTH ratio over the run separates deletion (reduced,
        unique seq) / segdup-masked deletion (normal + inside segdup) / copy-neutral
        LOH (normal + outside segdup); and
      * het BACKGROUND -- the flanks must be het-rich AND the run must be longer
        than expected under the local het rate (geometric tail < alpha), else it
        is a long ROH block, not a focal deletion.
    """
    positions = [p for (p, _r, _x) in lst]
    depths = [r + x for (_p, r, x) in lst]
    base = local_baselines(positions, depths, a.baseline_window_bp)

    # Per-site: het (anchors background) vs hom (homozygous-looking / hemizygous),
    # plus the local depth ratio. Zygosity is judged depth-INDEPENDENTLY.
    states = []
    for (pos, ref, alt), b in zip(lst, base):
        dp = ref + alt
        baf = alt / dp if dp else 0.0
        dr = dp / b if b else 1.0
        st = "het" if a.het_lo <= baf <= a.het_hi else "hom"
        states.append((pos, st, dr))
    n = len(states)

    def flank_het_rate(lo_idx, hi_idx):
        """het fraction among sites within +/- flank_bp outside [lo_idx, hi_idx]."""
        lp, rp = states[lo_idx][0], states[hi_idx][0]
        win = []
        k = lo_idx - 1
        while k >= 0 and lp - states[k][0] <= a.flank_bp:
            win.append(states[k]); k -= 1
        k = hi_idx + 1
        while k < n and states[k][0] - rp <= a.flank_bp:
            win.append(states[k]); k += 1
        if not win:
            return 0.0, 0, 0
        het = sum(1 for _p, st, _d in win if st == "het")
        return het / len(win), het, len(win)

    i = 0
    while i < n:
        if states[i][1] != "hom":
            i += 1
            continue
        # Extend a homozygous (LOH) run, tolerating up to max_het_in_run interior
        # het sites (a segdup paralog can leave a residual het inside a real del).
        j = i
        n_hom = het_in = 0
        last_hom = i
        drs = []
        while j < n:
            st = states[j][1]
            if st == "hom":
                n_hom += 1
                last_hom = j
                drs.append(states[j][2])
                j += 1
            elif het_in < a.max_het_in_run:
                het_in += 1
                j += 1
            else:
                break
        if n_hom >= a.min_hom_sites:
            start_pos, end_pos = states[i][0], states[last_hom][0]
            span = end_pos - start_pos
            mean_dr = sum(drs) / len(drs)
            fhr, fhet, fn = flank_het_rate(i, last_hom)
            cn = gcnv_cn_over(gcnv, chrom, start_pos, end_pos)

            # Is the run surprisingly long given the local het background?
            # geometric tail: P(>= n_hom homs in a row) ~ (1 - fhr)^n_hom
            p_loh = (1.0 - fhr) ** n_hom if fhr > 0 else 1.0
            het_flanked = (fn >= a.min_flank_sites and fhet >= a.min_flank_het
                           and fhr >= a.min_flank_het_rate and p_loh <= a.loh_alpha)
            focal = span < a.roh_min_bp
            depth_reduced = mean_dr < a.del_depth_max
            in_segdup = overlaps_bed(segdup, chrom, start_pos, end_pos)

            if not het_flanked:
                cls = "ROH_or_homtract"          # no het background -> not focal del
            elif not focal:
                cls = "large_LOH"                # Mb-scale -> not a focal rescue
            elif depth_reduced and cn is not None and cn < 2:
                cls = "gcnv_concordant_del"      # depth caller already saw it
            elif depth_reduced:
                cls = "depth_silent_del_candidate"   # unique-seq del depth missed
            elif in_segdup:
                cls = "segdup_masked_del_candidate"  # <-- the alpha-globin rescue
            else:
                cls = "CN_LOH_candidate"         # hemizygous, full depth, not segdup

            if cls in ("depth_silent_del_candidate", "gcnv_concordant_del",
                       "segdup_masked_del_candidate"):
                conf = ("high" if (n_hom >= a.high_hom_sites
                                   and fhr >= a.high_flank_het_rate
                                   and (mean_dr < a.high_depth_max or in_segdup))
                        else "medium" if n_hom >= a.min_hom_sites else "low")
            else:
                conf = "low"

            yield (chrom, start_pos, end_pos, span, n_hom, round(mean_dr, 3),
                   fhet, round(fhr, 3), round(p_loh, 4),
                   "." if cn is None else cn, cls, conf)
        i = max(j, i + 1)


def write_vcf(rows, sample, path):
    """Emit candidate rescues as a minimal annotation VCF (flag, not a call)."""
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##ALT=<ID=DEL,Description="Deletion">\n')
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        out.write('##INFO=<ID=HOMSITES,Number=1,Type=Integer,Description="consecutive homozygous common-SNP sites in the LOH run">\n')
        out.write('##INFO=<ID=DEPTHRATIO,Number=1,Type=Float,Description="mean local unique-read depth ratio over the run">\n')
        out.write('##INFO=<ID=FLANKHETRATE,Number=1,Type=Float,Description="het fraction in the flanking common-SNP windows (LOH-gate)">\n')
        out.write('##INFO=<ID=LOHP,Number=1,Type=Float,Description="geometric-tail p-value of the LOH run given the flank het rate">\n')
        out.write('##INFO=<ID=GCNVCN,Number=1,Type=String,Description="min gCNV copy number over the run (. if uncalled)">\n')
        out.write('##INFO=<ID=RESCUECLASS,Number=1,Type=String,Description="allelic rescue class">\n')
        out.write('##INFO=<ID=RESCUECONF,Number=1,Type=String,Description="rescue confidence">\n')
        out.write('##FILTER=<ID=AllelicRescue,Description="Focal hemizygous depth-reduced common-SNP run in a het background: deletion candidate the depth caller did not see (segdup-masked). Annotation/flag for review, not an auto-call.">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for (chrom, s, e, span, ns, dr, _fhet, fhr, ploh, cn, cls, conf) in rows:
            info = ("END=%d;SVTYPE=DEL;SVLEN=-%d;HOMSITES=%d;DEPTHRATIO=%.3f;"
                    "FLANKHETRATE=%.3f;LOHP=%.4f;GCNVCN=%s;RESCUECLASS=%s;RESCUECONF=%s"
                    % (e, span, ns, dr, fhr, ploh, cn, cls, conf))
            out.write("%s\t%d\t.\tN\t<DEL>\t.\tAllelicRescue\t%s\n"
                      % (chrom, s, info))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--allelic-counts", required=True,
                    help="GATK CollectAllelicCounts TSV (common SNPs, MQ>=30, recal CRAM)")
    ap.add_argument("--gcnv-vcf", default="",
                    help="gCNV genotyped_segments VCF for concordance annotation (optional)")
    ap.add_argument("--segdup-bed", default="",
                    help="segmental-duplication BED: focal normal-depth LOH inside a "
                         "segdup is classed segdup_masked_del_candidate (optional)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-vcf", default="",
                    help="optional candidate-DEL annotation VCF")
    # per-site
    ap.add_argument("--min-dp", type=int, default=8,
                    help="drop sites below this raw depth (noise floor)")
    ap.add_argument("--het-lo", type=float, default=0.25)
    ap.add_argument("--het-hi", type=float, default=0.75)
    # local baseline (refinement #1)
    ap.add_argument("--baseline-window-bp", type=int, default=100_000,
                    help="tile width for the local rolling-median depth baseline")
    ap.add_argument("--del-depth-max", type=float, default=0.75,
                    help="mean local depth ratio below this gates the LOH run as a deletion")
    # LOH-run building (PRIMARY signal)
    ap.add_argument("--min-hom-sites", type=int, default=5,
                    help="consecutive homozygous common SNPs needed to call an LOH run")
    ap.add_argument("--max-het-in-run", type=int, default=1,
                    help="interior het sites tolerated inside an LOH run (segdup residual)")
    ap.add_argument("--roh-min-bp", type=int, default=1_000_000,
                    help="LOH runs >= this span are reported as large_LOH, not focal dels")
    # LOH-gating against ROH (refinement #2)
    ap.add_argument("--flank-bp", type=int, default=50_000,
                    help="bp each side used to measure the het background")
    ap.add_argument("--min-flank-sites", type=int, default=4,
                    help="common SNPs required in the flank windows to judge background")
    ap.add_argument("--min-flank-het", type=int, default=2,
                    help="het sites required in the flanks (else ROH-like, gated out)")
    ap.add_argument("--min-flank-het-rate", type=float, default=0.15,
                    help="het fraction required in the flanks (LOH-gate vs ROH)")
    ap.add_argument("--loh-alpha", type=float, default=0.05,
                    help="geometric-tail p-value: run must be surprising given flank het rate")
    # confidence
    ap.add_argument("--high-hom-sites", type=int, default=7)
    ap.add_argument("--high-depth-max", type=float, default=0.60)
    ap.add_argument("--high-flank-het-rate", type=float, default=0.25)
    a = ap.parse_args()

    sites = parse_allelic_counts(a.allelic_counts, a.min_dp)
    gcnv = parse_gcnv_segments(a.gcnv_vcf)
    segdup = load_bed(a.segdup_bed)

    rows = []
    for chrom, lst in sites.items():
        rows.extend(scan_chrom(chrom, lst, a, gcnv, segdup))
    rows.sort(key=lambda r: (r[0], r[1]))

    with open(a.out_tsv, "w") as out:
        out.write("sample\tchrom\tstart\tend\tspan_bp\thom_sites\tmean_depth_ratio"
                  "\tflank_het\tflank_het_rate\tloh_p\tgcnv_CN\tclass\tconfidence\n")
        for r in rows:
            out.write(a.sample + "\t" + "\t".join(str(x) for x in r) + "\n")

    DEL_CLASSES = ("depth_silent_del_candidate", "gcnv_concordant_del",
                   "segdup_masked_del_candidate")
    if a.out_vcf:
        write_vcf([r for r in rows if r[10] in DEL_CLASSES], a.sample, a.out_vcf)

    n_resc = sum(1 for r in rows if r[10] == "depth_silent_del_candidate")
    n_segd = sum(1 for r in rows if r[10] == "segdup_masked_del_candidate")
    n_conc = sum(1 for r in rows if r[10] == "gcnv_concordant_del")
    n_cnloh = sum(1 for r in rows if r[10] == "CN_LOH_candidate")
    n_gate = sum(1 for r in rows if r[10] in ("ROH_or_homtract", "large_LOH"))
    print("%s: %d LOH runs (%d depth-silent del, %d segdup-masked del, "
          "%d gCNV-concordant del, %d CN-LOH, %d ROH/large gated)"
          % (a.sample, len(rows), n_resc, n_segd, n_conc, n_cnloh, n_gate))


if __name__ == "__main__":
    main()
