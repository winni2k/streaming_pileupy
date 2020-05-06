#!/usr/bin/env bash
# This script describes how the test fixtures in this file were created

# using samtools 1.9
mkdir -p fixture_1
pushd fixture_1
wd=$PWD
tmp=$(mktemp -d)
pushd $tmp
samtools split $wd/../base_1.sam
samtools mpileup -ABQ0 -f $wd/../base_1.fa *.bam > $wd/output.pileup
popd
rm -rf $tmp
samtools view -H $wd/../base_1.sam \
  | grep '^@RG' \
  | perl -pne 's/.*SM:(\S+).*/$1/' \
  > input.samples.txt
popd

mkdir -p fixture_2
pushd fixture_2
wd=$PWD
tmp=$(mktemp -d)
pushd $tmp
samtools split $wd/../base_1.sam
samtools mpileup -ABQ30 -f $wd/../base_1.fa *.bam > $wd/output.pileup
popd
rm -rf $tmp
samtools view -H $wd/../base_1.sam \
  | grep '^@RG' \
  | perl -pne 's/.*SM:(\S+).*/$1/' \
  > input.samples.txt
popd
