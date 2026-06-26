#!/bin/bash
# Breakpoint evidence for the HBA alpha-globin deletions (hg38), two complementary signals:
#
#  (A) SPAN-BASED junction reads. A split read whose two pieces (own alignment + its
#      SA supplementary) span a distance in the subtype's DELETION-SIZE range AND sit
#      at the right locus -> confirms + subtypes. Keying on the junction SPAN
#      (3.7 / 4.2 / 20 kb) is more robust than exact breakpoint coordinates, which
#      drift a few hundred bp between samples (microhomology, aligner choice).
#
#  (B) SOFT-CLIP PILEUP. Reads whose soft-clip BOUNDARY stacks within +/-W bp of a
#      canonical breakpoint, even with no SA tag. A pile of clips all ending at the
#      same base is breakpoint evidence -> sensitivity at weak junctions (low cov /
#      SA not emitted). Reported per 5'/3' breakpoint; a real junction piles at BOTH.
#
# Canonical hg38 breakpoints: -a3.7 173650<->177400 ; -a4.2 170000<->174200
# (segdup-internal) ; --SEA 165500<->185500. Usage: hba_breakpoints.sh <cram> <ref> <out.tsv>
set -uo pipefail
CRAM="$1"; REF="$2"; OUT="$3"
samtools view -T "$REF" "$CRAM" chr16:163000-187000 2>/dev/null | perl -e '
# subtype -> [5prime_bp, 3prime_bp, span_lo, span_hi]
my %sub=(
  "a3.7"=>[173650,177400,3400,4100],
  "a4.2"=>[170000,174200,3850,4600],
  "SEA" =>[165500,185500,17000,23000],
);
my $POSTOL=800;   # locus tolerance for span-based anchoring
my $W=60;         # clip-pileup tolerance (bp) around a canonical breakpoint
my $MINCLIP=15;   # minimum soft-clip length to count as a boundary
my %jc=("a3.7"=>0,"a4.2"=>0,"SEA"=>0,"other"=>0);
my %clip;         # binned clip-boundary position -> count
while(<STDIN>){
  my @f=split/\t/; next if @f<6;
  my ($flag,$pos,$cig)=($f[1],$f[3],$f[5]);
  next if ($flag & 0x100) || ($flag & 0x400);   # skip secondary / duplicate
  # (B) soft-clip boundary positions (10bp-binned)
  if($cig=~/^(\d+)S/ && $1>=$MINCLIP){ $clip{int($pos/10)*10}++ }
  if($cig=~/(\d+)S$/ && $1>=$MINCLIP){
    my $ref=0; while($cig=~/(\d+)([MDN=X])/g){ $ref+=$1 }
    $clip{int(($pos+$ref)/10)*10}++;
  }
  # (A) span-based SA junction
  next unless /SA:Z:chr16,(\d+),/;
  my $sa=$1; my $d=abs($pos-$sa);
  next if $d<2500 || $d>40000;
  my ($lo,$hi)=($pos<$sa?($pos,$sa):($sa,$pos));
  my $hit=0;
  for my $k (sort keys %sub){
    my($b5,$b3,$slo,$shi)=@{$sub{$k}};
    if($d>=$slo && $d<=$shi && abs($lo-$b5)<=$POSTOL && abs($hi-$b3)<=$POSTOL){ $jc{$k}++; $hit=1; last }
  }
  $jc{other}++ unless $hit;
}
sub pileup { my($bp)=@_; my $s=0; for my $p (keys %clip){ $s+=$clip{$p} if abs($p-$bp)<=$W } return $s }
for my $k ("a3.7","a4.2","SEA"){
  my($b5,$b3)=@{$sub{$k}};
  printf "%s_junction\t%d\n",$k,int($jc{$k}/2);   # /2: junction seen from both mates
  printf "%s_clip5\t%d\n",$k,pileup($b5);
  printf "%s_clip3\t%d\n",$k,pileup($b3);
}
printf "other_junction\t%d\n",int($jc{other}/2);
' > "$OUT"
