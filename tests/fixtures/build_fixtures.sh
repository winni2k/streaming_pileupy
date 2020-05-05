#!/usr/bin/env bash
# This script describes how the test fixtures in this file were created

# using samtools 1.9
cd fixture_1
samtools split fixture_1.sam
samtools mpileup -ABQ0 -h -OSAM *.bam > fixture_1.pileup
samtools view -H fixture_1.sam \
  | grep '^@RG' \
  | perl -pne 's/.*SM:(\S+).*/$1/' \
  | sort | uniq > fixture_1.samples.txt

