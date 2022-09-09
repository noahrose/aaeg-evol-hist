#!/bin/bash

cat sim.chr?.$1.$2*.viterbi | grep '1,1' | bedtools sort > sim.tracts.$1.$2.bed
cat sim2.chr?.$1.$2*.viterbi | grep '1,1' | bedtools sort > sim2.tracts.$1.$2.bed
bedtools coverage -a /tigress/noahr/ref/AaegL5_1mb_intervals.bed -b sim.tracts.$1.$2.bed -sorted > sim.tracts.$1.$2.1mb.bed
bedtools coverage -a /tigress/noahr/ref/AaegL5_1mb_intervals.bed -b sim2.tracts.$1.$2.bed -sorted > sim2.tracts.$1.$2.1mb.bed
