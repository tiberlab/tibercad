#!/bin/bash

if test $# -ne 1
then
  echo "Need output directory as argument" >&2
  exit 2
fi

dir=$1
mkdir -p $dir

echo "# angle P exx eyy ezz exz" > ${dir}/P_vs_angle.dat

for b in {0..90}
do
  bb=$(echo "$b+90" | bc -l)

  sed -e "s/@BETA@/$bb/" -e "s/@OUT@/$dir/" input.tib.in > input.tib
  ../../bin/tibercad input.tib

  p=$(awk '{ if (NR == 220) {P = $8;} if (NR == 250) { P = $8-P; print P; } }' ${dir}/dd_msh.dat)
  e="$(awk '{ if (NR == 250) { print $2, $3, $4, $7; } }' ${dir}/strain_msh.dat)"

  echo "$b $p $e" >> ${dir}/P_vs_angle.dat

done
