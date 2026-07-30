#!/usr/bin/env python3
"""HBA (alpha-globin) classifier for the GATK-gCNV postprocess pipeline.

Reports the three HBA copy-number compartments in the same schema as the Sentieon
segdup-caller (a3.7 / a4.2 / non-duplication) plus the alpha-thalassemia clinical
interpretation, and emits a YAML block compatible with that tool's output. The
interpretation logic is ported verbatim from segdup-caller genecaller/genes/hba.py.

Key design fact (validated empirically against the Sentieon reference): DEPTH-based
gCNV is blind to the paralog-masked heterozygous single-gene deletions (-a3.7/-a4.2)
because the ~95%%-identical surviving paralog back-fills the depth dip. So this
classifier FUSES THREE INDEPENDENT AXES rather than trusting depth alone:

  1. DEPTH-CN   -- normalized read depth over each region (catches large deletions,
                   the non-duplication backbone, and triplications; reliable in
                   unique sequence, blind in the paralog core);
  2. ALLELIC-LOH -- a focal run of hemizygous common SNPs over the locus (catches the
                   depth-silent het deletions that gCNV misses); and
  3. DOSAGE-DIRECTION -- a3.7-window vs a4.2-window relative depth, to lean the
                   subtype (HBA2-side => -a4.2, HBA1/3.7-side => -a3.7).

The integrated CN per gene compartment is depth-CN when depth is informative, else
an allelic-inferred call (marked, low-confidence) when depth reads 2 but the allelic
axis shows hemizygosity. Sentieon (when its YAML is supplied) is treated as ONE
REFERENCE among the axes -- NOT ground truth -- and every disagreement is reported.
"""
import argparse
import os
import re

# HBA region coordinates (hg38), from the segdup-caller gene model.
# a3.7 and a4.2 overlap (172871-174075); for dosage DIRECTION we use the
# NON-overlapping diagnostic sub-regions so the deleted gene is distinguishable.
REGIONS = {
    "a3.7": (172871, 176674),       # -a3.7 footprint
    "a4.2": (169818, 174075),       # -a4.2 footprint
    "a3.7_spec": (174075, 176674),  # -a3.7-only (HBA1 side) -> direction
    "a4.2_spec": (169818, 172871),  # -a4.2-only (HBA2 side) -> direction
    "locus": (169454, 177522),
}
REGION_COORDS = {
    "a3.7": "chr16:172871-176674",
    "a4.2": "chr16:169818-174075",
    "non-duplication": "Non-duplication region",
}
# dp_norm: cohort-median of (region_depth / genome_baseline) in a normal diploid
# (the value that should read CN=2). Calibrated on a 31-sample WGS cohort. non-dup
# is measured on the UNIQUE flanks (not the gene footprints) so a focal single-gene
# deletion does not spuriously drag the backbone CN down.
DP_NORM = {"a3.7_spec": 1.155, "a4.2_spec": 0.908, "nondup_flank": 1.005}
BP_CONFIRM_SPLIT = 8                # split-read count that confirms a breakpoint
HET_LO, HET_HI = 0.25, 0.75

# --------------------------------------------------------------- PSV dosage
# Paralog-specific-variant (PSV) HBA2/HBA1 dosage axis. At positions where the
# two duplication units differ in the reference, pool read support for the
# HBA2-defining base vs the HBA1-defining base across BOTH loci -> the pooled
# ratio estimates HBA2:HBA1 copy dosage, immune to the paralog cross-mapping
# that makes -a4.2 DEPTH-silent (the HBA2 gene body is ~4-copy near-identical
# sequence, so a single-copy loss barely dents depth). NOTE: this axis is
# advisory only -- it is NOT call-changing. It is weak for -a4.2 specifically
# (interlocus gene conversion homogenises the PSVs; few sites) so it must never
# be labelled "-a4.2": it is a generic HBA2-dosage-anomaly flag. See the cohort
# validation: it does NOT separate the four reference -a4.2 carriers, but it
# does surface real HBA2-depletion outliers every other axis calls normal.
# (hba2_pos, hba1_pos, hba2_base, hba1_base) -- 12 calibrated, isolated PSVs.
PSV_TABLE = [
    (173384, 177188, 'T', 'G'), (173615, 177426, 'A', 'G'),
    (173619, 177430, 'G', 'A'), (173621, 177432, 'T', 'G'),
    (173623, 177434, 'C', 'T'), (173626, 177437, 'C', 'T'),
    (173647, 177458, 'G', 'A'), (173649, 177460, 'G', 'C'),
    (173662, 177473, 'C', 'T'), (173664, 177475, 'T', 'C'),
    (173679, 177491, 'C', 'G'), (173726, 177538, 'G', 'A'),
]
# Sentieon-INDEPENDENT robust baseline: median + 1.4826*MAD of the HBA2 fraction
# over a 38-sample WGS cohort (label-free, same calibration style as DP_NORM).
PSV_BASELINE_MEDIAN = 0.522
PSV_BASELINE_SIGMA = 0.050
PSV_FLAG_Z = -2.0                   # robust-z below which an HBA2-dosage anomaly is flagged
PSV_LOWCOV_READS = 500             # pooled PSV reads below which the flag is low-confidence

_PILEUP_INDEL = re.compile(r'[+-](\d+)')
def _pileup_bases(ref_base, bases):
    """Resolve an mpileup base string to a list of called bases (. , -> ref)."""
    out, i, s = [], 0, bases
    while i < len(s):
        c = s[i]
        if c == '^':                # read-start marker + mapq char
            i += 2; continue
        if c in '$*':
            i += 1; continue
        if c in '+-':               # indel: skip the inserted/deleted run
            m = _PILEUP_INDEL.match(s[i:]); L = int(m.group(1))
            i += 1 + len(m.group(1)) + L; continue
        if c in '.,':
            out.append(ref_base); i += 1; continue
        if c in 'ACGTNacgtn':
            out.append(c.upper()); i += 1; continue
        i += 1
    return out

def psv_dosage(pileup_path):
    """Pool HBA2-base vs HBA1-base reads over PSV_TABLE. -> (hba2_fraction, total_reads)."""
    if not pileup_path or not os.path.exists(pileup_path):
        return None, 0
    want = set()
    for hp, lp, _, _ in PSV_TABLE:
        want.add(hp); want.add(lp)
    pile = {}
    for ln in open(pileup_path):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 5 and f[1].isdigit() and int(f[1]) in want:
            pile[int(f[1])] = _pileup_bases(f[2].upper(), f[4])
    h2 = h1 = 0
    for hp, lp, hb, lb in PSV_TABLE:
        for pos in (hp, lp):
            for b in pile.get(pos, []):
                if b == hb:
                    h2 += 1
                elif b == lb:
                    h1 += 1
    tot = h2 + h1
    return (h2 / tot if tot else None), tot


# --------------------------------------------------------------- depth -> CN
def depth_cn(region_dp, baseline_dp, dp_norm, baseline_cn=2):
    if baseline_dp <= 0:
        return None
    raw = baseline_cn * (region_dp / baseline_dp) / dp_norm
    return max(0, min(6, round(raw))), round(raw, 3)


def read_depths(path):
    """region-depth TSV: 'region<TAB>mean_depth' (incl. 'base')."""
    d = {}
    for ln in open(path):
        p = ln.split()
        if len(p) >= 2:
            try:
                d[p[0]] = float(p[1])
            except ValueError:
                pass
    return d


# --------------------------------------------------------------- allelic LOH
def allelic_evidence(counts_path, min_dp=8):
    """Per-region het/hom tally over the HBA locus from CollectAllelicCounts TSV.

    Returns dict region -> (n_sites, n_het, hemizygous_fraction, mean_depth).
    A hemizygous deletion drives het->hom (loss of het) over its footprint.
    """
    sites = []
    if counts_path and os.path.exists(counts_path):
        for ln in open(counts_path):
            if ln.startswith("@") or ln.startswith("CONTIG"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] != "chr16":
                continue
            pos, ref, alt = int(f[1]), int(f[2]), int(f[3])
            dp = ref + alt
            if dp >= min_dp and 169000 <= pos <= 178000:
                sites.append((pos, alt / dp, dp))
    out = {}
    for name, (s, e) in REGIONS.items():
        rs = [(baf, dp) for (p, baf, dp) in sites if s <= p <= e]
        n = len(rs)
        het = sum(1 for baf, _ in rs if HET_LO <= baf <= HET_HI)
        hemi = 1 - het / n if n else 0.0
        mdp = sum(dp for _, dp in rs) / n if n else 0.0
        out[name] = (n, het, round(hemi, 3), round(mdp, 1))
    return out


# --------------------------------------------------- clinical interpretation
# Ported verbatim from segdup-caller genecaller/genes/hba.py _get_clinical_interpretation
def clinical_interpretation(cn_a37, cn_a42, cn_nondup):
    if cn_nondup == 0:
        if cn_a37 == 0 and cn_a42 == 0:
            return "Homozygous large deletion (--/--) - Hemoglobin Bart's hydrops fetalis (lethal)"
        return "Abnormal non-duplication CN=0 - Requires manual review"
    if cn_nondup == 1:
        if cn_a37 == 0 and cn_a42 == 0:
            return "Heterozygous large deletion, likely --SEA, --MED, --FIL, or --THAI (αα/--) - Silent carrier"
        elif cn_a37 == 1 and cn_a42 == 2:
            return "Compound -α3.7 with large deletion (-α3.7/--) - Hemoglobin H disease"
        elif cn_a37 == 2 and cn_a42 == 1:
            return "Compound -α4.2 with large deletion (-α4.2/--) - Hemoglobin H disease"
        return "Abnormal non-duplication CN=1 with unusual pattern - Requires manual review"
    if cn_a37 == 2 and cn_a42 == 2:
        return "Normal/Reference (αα/αα) - Normal"
    if cn_a37 == 1 and cn_a42 == 2:
        return "Heterozygous -α3.7 deletion (αα/-α3.7) - Silent carrier"
    if cn_a37 == 2 and cn_a42 == 1:
        return "Heterozygous -α4.2 deletion (αα/-α4.2) - Silent carrier"
    if cn_a37 == 0 and cn_a42 == 2:
        return "Homozygous -α3.7 deletion (-α3.7/-α3.7) - Alpha-thalassemia trait"
    if cn_a37 == 2 and cn_a42 == 0:
        return "Homozygous -α4.2 deletion (-α4.2/-α4.2) - Alpha-thalassemia trait"
    if cn_a37 == 1 and cn_a42 == 1:
        return "Compound heterozygous (-α3.7/-α4.2) - Alpha-thalassemia trait"
    if cn_a37 == 3 and cn_a42 == 2:
        return "α-globin triplication involving -α3.7 region - Usually normal"
    if cn_a37 == 2 and cn_a42 == 3:
        return "α-globin triplication involving -α4.2 region - Usually normal"
    if cn_nondup > 2:
        return f"Complex rearrangement (non-dup CN={cn_nondup}) - Requires manual review"
    return (f"Unusual copy number pattern (-α3.7={cn_a37}, -α4.2={cn_a42}, "
            f"non-dup={cn_nondup}) - Requires manual review")


DISPLAY = {"a3.7": "-α3.7", "a4.2": "-α4.2", "non-duplication": "non-duplication"}


def region_line(name, cn, inferred=False):
    disp = DISPLAY.get(name, name)
    coords = REGION_COORDS.get(name, REGION_COORDS["non-duplication"])
    tail = " [allelic-inferred, depth-silent]" if inferred else ""
    if cn == 0:
        return f"{disp}: Homozygous deletion (CN=0, {coords}){tail}"
    if cn == 1:
        return f"{disp}: Heterozygous deletion (CN=1, {coords}){tail}"
    if cn == 2:
        return f"{disp}: Normal copy number (CN=2)"
    if cn == 3:
        return f"{disp}: Triplication (CN=3, {coords}){tail}"
    return f"{disp}: Copy number variation (CN={cn}, {coords}){tail}"


# ------------------------------------------------------- Sentieon reference
def read_sentieon_hba(path):
    """Extract a3.7/a4.2/non-duplication CN from a Sentieon segdup-caller YAML."""
    if not path or not os.path.exists(path):
        return None
    cns, in_hba, in_cn = {}, False, False
    for ln in open(path):
        if ln.startswith("HBA:"):
            in_hba = True
            continue
        if in_hba and ln and not ln[0].isspace():
            break
        if in_hba and ln.strip() == "Copy numbers:":
            in_cn = True
            continue
        if in_cn:
            if ln.strip().endswith(":") or (ln.strip() and not ln.startswith("    ")):
                in_cn = False
            else:
                p = ln.strip().split(":")
                if len(p) == 2 and p[1].strip().lstrip("-").isdigit():
                    cns[p[0].strip()] = int(p[1])
    return cns or None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--depths", required=True, help="region-depth TSV (region<TAB>mean_dp, incl 'base')")
    ap.add_argument("--allelic-counts", default="", help="CollectAllelicCounts TSV over the HBA locus")
    ap.add_argument("--breakpoints", default="", help="breakpoint-read stats TSV (split/bigclip/disc<TAB>count)")
    ap.add_argument("--psv-pileup", dest="psv_pileup", default="", help="samtools mpileup over the HBA2/HBA1 PSV loci (advisory HBA2-dosage axis)")
    ap.add_argument("--variants-vcf", default="",
                    help="provenance only: path recorded verbatim in the YAML 'Variants' "
                         "field. Not parsed, and not necessarily a VCF -- the caller "
                         "passes the per-sample allelic-counts TSV that backs the call.")
    ap.add_argument("--sentieon-yaml", default="", help="optional Sentieon YAML to cross-check (reference, not truth)")
    ap.add_argument("--out-yaml", required=True)
    ap.add_argument("--out-crosscheck", required=True)
    ap.add_argument("--hemi-call", type=float, default=0.80,
                    help="hemizygous-fraction over a region above which allelic axis calls LOH")
    ap.add_argument("--junc-min", type=int, default=2,
                    help="junction-spanning split reads needed to confirm a recurrent junction")
    ap.add_argument("--clip-min", type=int, default=4,
                    help="soft-clip pileup at BOTH breakpoints needed to confirm a junction (sensitivity)")
    a = ap.parse_args()

    dp = read_depths(a.depths)
    base = dp.get("base", 0.0)
    # JUNCTION-ANCHORED breakpoint evidence (paralog-immune; subtype-specific).
    # Counts of split reads spanning each recurrent alpha-globin junction, with the
    # SA mate required on the partner breakpoint -- rejects segdup-edge artifacts and
    # names the lesion by which junction is spanned.
    junc = {}
    if a.breakpoints and os.path.exists(a.breakpoints):
        for ln in open(a.breakpoints):
            p = ln.split()
            if len(p) >= 2:
                try:
                    junc[p[0]] = int(float(p[1]))
                except ValueError:
                    pass
    jget = lambda k: junc.get(k, 0)
    # A subtype junction is confirmed by EITHER >= junc_min span-based SA reads OR a
    # soft-clip PILEUP at BOTH its canonical breakpoints (sensitivity where the SA
    # tag is not emitted; requiring both ends keeps it specific).
    def confirmed(sub):
        return (jget(sub + "_junction") >= a.junc_min or
                (jget(sub + "_clip5") >= a.clip_min and jget(sub + "_clip3") >= a.clip_min))
    a37_jc, a42_jc, sea_jc = confirmed("a3.7"), confirmed("a4.2"), confirmed("SEA")
    any_jc = a37_jc or a42_jc or sea_jc or jget("other_junction") >= a.junc_min

    # depth-CN per diagnostic region (specific = non-overlapping, for direction)
    dcn = {}
    for name in ("a3.7_spec", "a4.2_spec"):
        dcn[name] = depth_cn(dp.get(name, 0.0), base, DP_NORM[name]) if name in dp else None
    nondup_flank = (dp.get("flank5", 0.0) + dp.get("flank3", 0.0)) / 2
    dcn_nd = depth_cn(nondup_flank, base, DP_NORM["nondup_flank"])
    allelic = allelic_evidence(a.allelic_counts)

    # dosage direction from the NON-overlapping sub-regions (the deleted gene is lower)
    raw37 = dcn["a3.7_spec"][1] if dcn["a3.7_spec"] else None
    raw42 = dcn["a4.2_spec"][1] if dcn["a4.2_spec"] else None
    lean = "unresolved"
    if raw37 is not None and raw42 is not None:
        if raw42 < raw37 - 0.12:
            lean = "-a4.2 (HBA2 lower)"
        elif raw37 < raw42 - 0.12:
            lean = "-a3.7 (HBA1 lower)"
        else:
            lean = "symmetric"

    # ---- integrate depth + breakpoint (conservative: subtype only from a depth
    #      dip; breakpoint flags an SV without forcing a noisy subtype guess) ----
    # NOTE: an allelic-inferred CN path (flip `inferred`->True when depth reads 2
    # but allelic shows hemizygosity) is intentionally NOT wired: the BAF/LOH axis
    # is saturated (~1.0) cohort-wide over the alpha segdup and is non-discriminating
    # for -a4.2 (validated). `--hemi-call`/`loc_hemi` are reported as evidence only.
    # The HBA2-dosage signal lives on the PSV axis below (advisory flag, not a call).
    inferred = {"a3.7": False, "a4.2": False}
    conf = {"a3.7": "n/a", "a4.2": "n/a", "non-duplication": "n/a"}
    loc_hemi = allelic.get("locus", (0, 0, 0.0, 0.0))[2]
    cn = {}
    cn["non-duplication"] = dcn_nd[0] if dcn_nd else 2
    conf["non-duplication"] = "depth(flank)" if dcn_nd else "n/a"
    jc_of = {"a3.7": a37_jc, "a4.2": a42_jc}
    for name, spec in (("a3.7", "a3.7_spec"), ("a4.2", "a4.2_spec")):
        d = dcn[spec][0] if dcn[spec] else 2
        raw = dcn[spec][1] if dcn[spec] else 2.0
        depth_loss = raw < 1.6           # clear dip (not just rounding)
        if d < 2 or (jc_of[name] and depth_loss):
            cn[name] = 1
            conf[name] = ("high (depth-dip+junction)" if (d < 2 and jc_of[name])
                          else "high (depth-dip)" if d < 2
                          else "high (junction+depth-loss)")
        else:
            cn[name] = 2
            conf[name] = "depth-normal" + (" (junction present, no depth loss -> not a deletion)"
                                           if jc_of[name] else "")
    # --SEA: a confirmed large-deletion junction sets the non-dup backbone to 1.
    if sea_jc:
        cn["non-duplication"] = min(cn["non-duplication"], 1)
        conf["non-duplication"] = "junction-confirmed (--SEA)"
    cn_a37, cn_a42, cn_nd = cn["a3.7"], cn["a4.2"], cn["non-duplication"]
    # -a4.2 is depth/junction/BAF/PSV-silent (HBA2 footprint is ~4-copy near-identical
    # paralog sequence). A CN=2 here is NOT a negative -- it means "not assessable".
    if cn_a42 == 2:
        conf["a4.2"] = ("not assessable by WGS (HBA2 footprint ~4-copy paralog; "
                        "depth/junction/BAF/PSV-silent) - -a4.2 NOT excluded")

    # A junction-confirmed SV that is copy-number-neutral on depth is a REAL
    # rearrangement (artifact-filtered) below depth's subtype resolution -- e.g. a
    # tandem-duplication junction (triplication, gain) or a balanced/complex event.
    sv_review = any_jc and cn_a37 == 2 and cn_a42 == 2 and cn_nd == 2
    sv_kind = ("tandem-duplication junction (triplication-suspect)"
               if (jget("a4.2_junction") >= a.junc_min or jget("other_junction") >= a.junc_min)
               else "junction-confirmed rearrangement")

    # ---- advisory PSV HBA2-dosage axis (NOT call-changing) ----
    psv_frac, psv_reads = psv_dosage(a.psv_pileup)
    psv_z = (psv_frac - PSV_BASELINE_MEDIAN) / PSV_BASELINE_SIGMA if psv_frac is not None else None
    psv_flag = psv_z is not None and psv_z < PSV_FLAG_Z
    psv_conf = "low(thin-coverage)" if psv_reads < PSV_LOWCOV_READS else "ok"

    # ---- interpretation (ported schema) ----
    lines = [region_line("non-duplication", cn_nd),
             region_line("a3.7", cn_a37, inferred["a3.7"]),
             region_line("a4.2", cn_a42, inferred["a4.2"])]
    clin = clinical_interpretation(cn_a37, cn_a42, cn_nd)
    if sv_review:
        clin = (f"Structural variant at HBA - {sv_kind} (junction-confirmed: "
                f"a3.7={jget('a3.7_junction')}, a4.2={jget('a4.2_junction')}, "
                f"SEA={jget('SEA_junction')}, other={jget('other_junction')}) but "
                f"copy-number-neutral on depth - Requires manual review")
    # -a4.2 is never assessable on WGS: a CN=2 does NOT exclude it.
    if cn_a42 == 2:
        clin += (" | NOTE: -a4.2 not assessable by WGS (paralog-masked) - NOT excluded; "
                 "orthogonal assay (gap-PCR/MLPA) required to rule out")
    # Possible compound genotype: a depth-visible -a3.7 plus an independent a4.2
    # junction or HBA2-dosage depletion (the -a3.7/-a4.2 compound pattern, e.g. ZHAO-CY).
    if cn_a37 == 1 and (a42_jc or psv_flag):
        clin += (" | REVIEW: possible second HBA lesion (compound -a3.7/-a4.2?) - "
                 f"{'a4.2 junction confirmed' if a42_jc else ''}"
                 f"{' + ' if (a42_jc and psv_flag) else ''}"
                 f"{('HBA2-dosage depleted (PSV robust_z=%.2f)' % psv_z) if psv_flag else ''}")
    # Generic HBA2-dosage anomaly (advisory; NOT labelled -a4.2).
    if psv_flag:
        clin += (f" | FLAG: HBA2-dosage anomaly (PSV robust_z={psv_z:.2f}, conf={psv_conf}) - "
                 "recommend orthogonal testing (gap-PCR/MLPA)")
    lines.append(f"Clinical interpretation: {clin}")
    interp = "\n".join(lines)

    # ---- emit YAML (HBA block, Sentieon-compatible + evidence) ----
    def y(s, ind=2):
        return " " * ind + s
    out = ["HBA:"]
    out.append(y("Copy numbers:"))
    out.append(y(f"a3.7: {cn_a37}", 4))
    out.append(y(f"a4.2: {cn_a42}", 4))
    out.append(y(f"non-duplication: {cn_nd}", 4))
    if a.variants_vcf:
        out.append(y(f"Variants: {a.variants_vcf}"))
    out.append(y("cn_interpretation: |-"))
    for l in interp.split("\n"):
        out.append(y(l, 4))
    out.append(y("evidence:"))
    out.append(y("caller: GATK-gCNV-postprocess + allelic-LOH + breakpoint (createpanelrefs HBA_CLASSIFIER)", 4))
    out.append(y("confidence:", 4))
    for k in ("a3.7", "a4.2", "non-duplication"):
        out.append(y(f"{k}: {conf[k]}", 6))
    out.append(y("depth_cn:", 4))
    for k in ("a3.7_spec", "a4.2_spec"):
        v = dcn[k]
        out.append(y(f"{k}: {v[0]} (raw {v[1]})" if v else f"{k}: null", 6))
    out.append(y(f"non-duplication_flank: {dcn_nd[0]} (raw {dcn_nd[1]})" if dcn_nd else "non-duplication_flank: null", 6))
    out.append(y(f"baseline_depth: {round(base,1)}", 4))
    out.append(y("allelic_loh:", 4))
    for k in ("a3.7", "a4.2", "locus"):
        n, het, hemi, mdp = allelic[k]
        out.append(y(f"{k}: {{sites: {n}, het: {het}, hemizygous_frac: {hemi}, mean_dp: {mdp}}}", 6))
    out.append(y(f"dosage_direction: {lean}", 4))
    out.append(y("breakpoint_junctions:", 4))
    for sub in ("a3.7", "a4.2", "SEA"):
        out.append(y(f"{sub}: {{span_reads: {jget(sub + '_junction')}, "
                     f"clip5: {jget(sub + '_clip5')}, clip3: {jget(sub + '_clip3')}, "
                     f"confirmed: {str(confirmed(sub)).lower()}}}", 6))
    out.append(y(f"other_junction: {jget('other_junction')}", 6))
    out.append(y(f"allelic_inferred: {{a3.7: {str(inferred['a3.7']).lower()}, a4.2: {str(inferred['a4.2']).lower()}}}", 4))
    if psv_frac is not None:
        out.append(y(f"psv_dosage: {{hba2_fraction: {round(psv_frac, 3)}, reads: {psv_reads}, "
                     f"cohort_median: {PSV_BASELINE_MEDIAN}, robust_z: {round(psv_z, 2)}, "
                     f"confidence: {psv_conf}}}", 4))
        if psv_flag:
            out.append(y(f'hba2_dosage_anomaly: "robust_z={round(psv_z, 2)} (conf={psv_conf}) '
                         '-> recommend orthogonal testing (gap-PCR/MLPA); advisory, NOT an -a4.2 call"', 4))
    else:
        out.append(y("psv_dosage: null  # no PSV pileup provided", 4))
    open(a.out_yaml, "w").write("\n".join(out) + "\n")

    # ---- multi-way cross-check (no axis is truth) ----
    sent = read_sentieon_hba(a.sentieon_yaml)
    cc = [f"# HBA cross-check  sample={a.sample}",
          f"# compartments: a3.7 / a4.2 / non-duplication",
          f"call(GATK_depth+allelic)\t{cn_a37}/{cn_a42}/{cn_nd}\t{clin}",
          f"axis_depth_only\t{dcn['a3.7_spec'][0] if dcn['a3.7_spec'] else '.'}/"
          f"{dcn['a4.2_spec'][0] if dcn['a4.2_spec'] else '.'}/{dcn_nd[0] if dcn_nd else '.'}",
          f"axis_allelic_locus_hemizygous_frac\t{loc_hemi}",
          f"axis_dosage_direction\t{lean}",
          f"axis_breakpoint_junction\ta3.7(span={jget('a3.7_junction')},clip={jget('a3.7_clip5')}/{jget('a3.7_clip3')}) "
          f"a4.2(span={jget('a4.2_junction')},clip={jget('a4.2_clip5')}/{jget('a4.2_clip3')}) "
          f"SEA(span={jget('SEA_junction')}) other={jget('other_junction')}"]
    if psv_frac is not None:
        cc.append(f"axis_psv_hba2_fraction\t{round(psv_frac, 3)} "
                  f"(robust_z={round(psv_z, 2)}, reads={psv_reads}, conf={psv_conf})")
        if psv_flag:
            cc.append(f"hba2_dosage_anomaly\trobust_z={round(psv_z, 2)} -> recommend orthogonal "
                      "testing (gap-PCR/MLPA) [advisory, NOT an -a4.2 call]")
    if sent is not None:
        sv = (sent.get("a3.7", "."), sent.get("a4.2", "."), sent.get("non-duplication", "."))
        agree = (sv[0] == cn_a37 and sv[1] == cn_a42 and sv[2] == cn_nd)
        cc.append(f"reference_sentieon\t{sv[0]}/{sv[1]}/{sv[2]}\t(reference, not truth)")
        cc.append(f"concordance_vs_sentieon\t{'CONCORDANT' if agree else 'DISCORDANT'}")
        if not agree:
            diffs = []
            for k, mine, theirs in (("a3.7", cn_a37, sv[0]), ("a4.2", cn_a42, sv[1]),
                                    ("non-dup", cn_nd, sv[2])):
                if mine != theirs:
                    diffs.append(f"{k}: mine={mine} sentieon={theirs}")
            cc.append("discordant_compartments\t" + "; ".join(diffs))
    else:
        cc.append("reference_sentieon\tNA")
    open(a.out_crosscheck, "w").write("\n".join(cc) + "\n")
    print(f"{a.sample}: HBA call {cn_a37}/{cn_a42}/{cn_nd} | {clin}"
          + (f" | vs Sentieon {sent.get('a3.7','.')}/{sent.get('a4.2','.')}/{sent.get('non-duplication','.')}"
             if sent else ""))


if __name__ == "__main__":
    main()
