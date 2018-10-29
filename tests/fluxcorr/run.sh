#!/usr/bin/env bash

tee input.par <<EOF
mbh 10
mdot 0.1
radius 10
alpha 0.1
eta 0.1
nu 0.25
EOF

mkdir -p data

zdisks="15 20 40 80 120 200 400 800 2000 4000 8000"
DVFLAGS='-compton -no-alpha-prad'
parallel cat input.par \| diskvert $DVFLAGS -htop {}    -fluxcorr -o data/f{} ::: $zdisks
parallel cat input.par \| diskvert $DVFLAGS -htop {} -no-fluxcorr -o data/z{} ::: $zdisks
parallel cat input.par \| diskvert $DVFLAGS    -fluxcorr -o data/fauto ::: $zdisks
parallel cat input.par \| diskvert $DVFLAGS -no-fluxcorr -o data/zauto ::: $zdisks

rm -f input.par

./plot.py $zdisks
