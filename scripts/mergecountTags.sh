#!/bin/bash
#
#$Author: Anthony Boureux <Anthony.boureux@univ-montp2.fr> $
#
# merge countTags values in one file

# default local variable
# default variable
TODAY=`date +%Y%m%d%H%M`

#Usage function to print help info.
usage () {
   cat <<EOF
   Concatenate all countTags files, and output do stdout

Usage: $0 "parameters"
with option
  -o file   : output in file instead of STDOUT
  -d dir	  : take all files from that directory, and from argument line
  -n        : kmer name are present in countTags file
  -s        : regroup summary files as well
  -v        : verbose mode
	-h				: help
EOF
   exit 1
}

# variables
# kmer name are present or not
countcol=2
# do you group summary file ?
dosummary=0
# output filename
outputfile=""

#Read in the various options.
while getopts o:d:nsvh OPT
do
	case $OPT in
    o)  outputfile=$OPTARG;;
		d)	dir=$OPTARG;;
    n)  countcol=3;;
    s)  dosummary=1;;
		v)	debug=1;;
		h)	usage;;
		\?)	echo "Wrong arguments"
			usage;;
	esac
done

shift `expr $OPTIND - 1`

# output to stderr
[ $debug ] && >&2 echo "dir = $dir"
[ $debug ] && >&2 echo "count column = $countcol"

# get list of files from $dir and $@
files=''

if [ $# ]
then
  files=$@
fi
if [ "$dir" != "" ]
then
  files="$files $(find $dir -type f -name '*.tsv')"
fi
[ $debug ] && >&2 echo "files= $files"

# create a tempory directory
temp=$(mktemp -d)

if [ ! -d $temp ]
then
  echo "Error can't create a tempory directory"
  echo "Check your permissions"
  exit 1
fi
[ $debug ] && >&2 echo "tmpdir= $temp"

# store filename of one file to get the first two columns
afile=''

# extract all column from $countcol
for file in $files
do
  cut -f $countcol- $file > $temp/$(basename $file)
  [ -s ${file/tsv/summary} ] && tail -n 3 ${file/tsv/summary} | cut -f $countcol-  > $temp/$(basename $file .tsv).summary
  afile=$file
done

[ $debug ] && >&2 echo "one file analysed = $afile"

# create the merged files
if [ ! -s $afile ]
then
  echo "No countTags found"
  exit 1
fi


if [ ! -s $temp/$(basename $afile) ]
then
  echo "No count column files found: $temp/$(basename $afile)"
  exit 1
fi
colname=$(($countcol - 1))
paste <(cut -f 1,$colname $afile) $temp/*.tsv > $outputfile.tsv
cat <(head -7 ${afile/tsv/summary}) <(paste <(tail -n 3 ${afile/tsv/summary} | cut -f 1,$colname) $temp/*.summary) > $outputfile.summary

if [ $? -eq 0 ]
then
  echo "Completed job" >&2
  rm -rf $temp
else
  echo "Error when merging all file from $temp"
  exit 1
fi
