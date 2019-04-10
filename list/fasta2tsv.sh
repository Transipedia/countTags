#!/bin/bash
# convert fasta to tsv format for countTags
# remove '>', replace ' ' in header by '_'
for f in $*
do
  paste <(grep -v '>' $f) <(grep  '>' $f | sed 's/>//; s/ /_/g') > ${f/fa/tsv}
done
