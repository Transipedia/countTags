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
  -o file   : output in file instead of STDOUT, without extension
  -d dir	  : take all files from that directory, and from argument line
  -n        : kmer name are present in countTags file
  -v        : verbose mode
	-h				: help
EOF
   exit 1
}

# variables
# kmer name are present or not
countcol=2
# output filename
outputfile=""

#Read in the various options.
while getopts o:d:nsvh OPT
do
	case $OPT in
    o)  outputfile=$OPTARG;;
		d)	dir=$OPTARG;;
    n)  countcol=3;;
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

# test if everthing was rigth
OKrm=0

# Do we output to STDOUT or in a file ?
outputformat="/dev/stdout"
if [ "$outputfile" != "" ]
then
  outputformat="$outputfile.tsv"
fi

paste <(cut -f 1,$colname $afile) $temp/*.tsv > $outputformat
if [ $? -eq 0 ]
then
  echo "Completed job for tsv files" >&2
  OKrm=1
else
  echo "Error when merging all tsv files from $temp"
  exit 1
fi

# Test if we have summary file
if [ -s ${afile/tsv/summary} ]
then
  OKrm=0   # turn off in case problem with summary files
  # take name from $dir if no $outputfile given
  if [ "$outputfile" == "" ]
  then
    if [ "$dir" != "" ]
    then
      outputfile=$(basename $dir)
    else
      outputfile=$TODAY-results
    fi
  fi
  cat <(head -7 ${afile/tsv/summary}) <(paste <(tail -n 3 ${afile/tsv/summary} | cut -f 1,$colname) $temp/*.summary) > $outputfile.summary
  if [ $? -eq 0 ]
  then
    echo "Completed job for summary files" >&2
    OKrm=1
  else
    echo "Error when merging all summary files from $temp"
    exit 1
  fi
fi

if [ "$OKrm" == 1 ]
then
  echo "Completed all jobs" >&2
  rm -rf $temp
else
  echo "Error when merging some files from $temp"
  exit 1
fi
