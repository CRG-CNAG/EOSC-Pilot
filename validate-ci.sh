#!/bin/bash 
#
# basic output files check 
# improve it for a stronger validation check
# 
[[ -s results/gonl-1a_L1.recal.sorted.bam ]] || false
[[ -s results/gonl-1a_L7.recal.sorted.bam ]] || false
