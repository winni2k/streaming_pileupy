#!/usr/bin/env bash
# This script describes how the test fixtures in this file were created

# using samtools 1.9
pushd fixture_1
wd=$PWD
tmp=$(mktemp -d)
pushd $tmp
samtools split $wd/fixture_1.sam
samtools mpileup -ABQ0 -f $wd/fixture_1.fa *.bam > $wd/fixture_1.pileup
popd
rm -rf $tmp
samtools view -H fixture_1.sam \
  | grep '^@RG' \
  | perl -pne 's/.*SM:(\S+).*/$1/' \
  > fixture_1.samples.txt

