# countTags - counting a set of k-mers from a set of FASTQ files

## Installation

1. Clone this repository : `git clone --recursive https://github.com/Transipedia/countTags.git`
    or git clone https://github.com/Transipedia/countTags.git && git submodule init && git submodule update
2. Compile the software : `make`
3. Place the binary in a directory which is in your `$PATH`

## Usage

>
> countTags -k 31 -i my_kmers.file.tsv data1.fastq.gz data2.fastq.gz ...
>

You need two files:
 * a file with the kmers that you want quantify in fastq data.
 * one or more fastq data file (compressed with gz or not)

### Input format for kmers

You create a file containing the tags/kmers you want to quantify.

This file can be in three tags/kmers differents formats:

 * fasta format
 * or separator format : sequence name (the separator can be a space or a tabulation or a comma or a semi-comma)
 * or raw format : only the tag/kmer sequence

* As fasta format:

```
>kmer.1
CACGTACTACGTTGTAGCCCACTTCCACTA
>kmer.2
GCGGGGTCGAAGAAGGTGGTGTTGAGGTTG
>kmer.3
GTTGGCCGAGTGGAGACTGGTGTTCTCAAA
>kmer.4
TGTTGCCATGGTAATCCTGCTCAGTACGAG
>kmer.5
GCTTAGGCAGAAGCCCTATTACTTTGCAAG
>kmer.6
ATAGGGGAAATCAGTGAATGAAGCCTCCTA
```

 * As csv/tsv format:

```
CACGTACTACGTTGTAGCCCACTTCCACTA;kmer.1
GCGGGGTCGAAGAAGGTGGTGTTGAGGTTG;kmer.2
GTTGGCCGAGTGGAGACTGGTGTTCTCAAA;kmer.3
TGTTGCCATGGTAATCCTGCTCAGTACGAG;kmer.4
GCTTAGGCAGAAGCCCTATTACTTTGCAAG;kmer.5
ATAGGGGAAATCAGTGAATGAAGCCTCCTA;kmer.6
```

 * As raw format:

```
CACGTACTACGTTGTAGCCCACTTCCACTA
GCGGGGTCGAAGAAGGTGGTGTTGAGGTTG
GTTGGCCGAGTGGAGACTGGTGTTCTCAAA
TGTTGCCATGGTAATCCTGCTCAGTACGAG
GCTTAGGCAGAAGCCCTATTACTTTGCAAG
ATAGGGGAAATCAGTGAATGAAGCCTCCTA
```

You must use option '-i' to specify this kmers file.  You can provide the kmer
file via the standard input by using '-i -' as filename.  If you file is
gziped, you can pass directly with the '-i mytags.gz' option or the pipe if
needed, but if it is in other compression format, uncompress the file with the
right tool and pass to countTags via the pipe and option '-i -'.

### kmers length

All tags/kmers must have at least the K-mer length, if too short, tags/kmers are discarded.
They are print to STDERR.

The maximum authorize tag length is 32 bp (one integer).

K-mer length can be provided to countTags using the `-k INT` option to change the default option = 31 (from version 0.6).

For example :

`countTags -k 31 file.fa file1.fastq.gz file2.fastq.gz`

### Managing stranded files

By default, countTags count the canonical tag/kmer between the forward and the reverse tag
(the first one in alphabetical order) and output this sequence in the result.

If you want only one strand to be count, you have to provide the tag sequence in the strand that
you want to count  and use the option '--stranded' for countTags.
Therefore, for stranded paired fastq files you will count only one pair, the one in the same
strand that your tag. For paired fastq see below the '--paired' option.

### Managing Paired-End files

You can now use the '--paired format' option to count stranded pair-end fastq file.
You have to specify the pair-end format: either 'rf' (the most used), 'fr' or 'rr'
in accordance with the library setup. In this case, only the two paired fastq file must be given.

When you set the paired option, the stranded option is set to true, otherwise this is
meaning nothing.

For paired-end files, you can use the '--merge-counts' option to get the total count for the sample.

### Normalize count values

For now, countTags can normalize the values of each tag/kmer with the option `-n|--normalize`.
In this case the values are millions of tag/kmer in each sample.

You can normalize by billions of tag/kmer using the option `-b|--billions`.
It will be the default normalization from version 1.0.

### Output reads

You can use option `-r|--reads` to output reads that are matched by a tag/kmer.

If more than one tag/kmer match a read, the read is output only once with all tag/kmer sequences and names separates by a comma.

The output format is a tabular file with :

 * all tag/kmer sequences matching the read
 * names of all tag/kmers if option `-t|--tag-names`
 * fastq file name
 * read header
 * read sequence
 * read quality line

A fastq file for each input fastq files are generated with the name:

 * name_given_with_option--reads + "-" + input fastq file.

```
  countTags -i test/TAGS_test.csv -k 30 --tag-names --reads /tmp/reads test/test_?.fastq.gz
  ls /tmp/reads*
  # /tmp/reads  /tmp/reads-test_1.fastq  /tmp/reads-test_2.fastq
```

## Bugs

 * 20210425-01 : Count twice the read if read are paired and overlapping
