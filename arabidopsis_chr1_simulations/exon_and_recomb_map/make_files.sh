#! /bin/bash

( cat sim_seq_info.txt | grep "recRate" ) > recomb.txt
( echo "start,stop"; cat sim_seq_info.txt | grep "exon" | awk 'BEGIN {OFS=","}{print $2,$3}') > exons.csv