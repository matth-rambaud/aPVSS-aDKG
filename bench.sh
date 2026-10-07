#!/bin/bash
nb=20

bench_pvss() {
for (( c=1; c<=nb; c++ ))
do
    echo
    echo -e "\e[32m#######  PVSS $1 $2 ############ ITERATION $c/$nb #################################\e[0m" 
    echo
    ./build/benchs/pvss_adaptative $1 $2
done
}
 
bench_pvss 50 24
bench_pvss 100 49
bench_pvss 150 74
bench_pvss 200 99


bench_dkg() {
for (( c=1; c<=nb; c++ ))
do
    echo
    echo -e "\e[32m########  DKG $1 $2 $3 ############ ITERATION $c/$nb #################################\e[0m" 
    echo

    ./build/benchs/dkg_adaptative $1 $2 $3
done
}
 
bench_dkg 10 4 0
bench_dkg 20 9 0
bench_dkg 30 14 0
bench_dkg 40 19 0
bench_dkg 50 24 0
bench_dkg 10 4 1
bench_dkg 20 9 1
bench_dkg 30 14 1
bench_dkg 40 19 1
bench_dkg 50 24 1

