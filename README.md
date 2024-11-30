# countTags - counting a set of k-mers into a set of FASTQ files

## Installation

1. Clone this repository : `git clone --recursive https://gitlab:get/counttags.git`
    or git clone https://gitlab:get/counttags.git && git submodule init && git submodule update
2. Compile the software : `make`
3. Place the binary in a directory which is in your `$PATH`

## Usage

First you need to create a file containing the tags/kmers you want to quantify.

You can use three tags/kmers format:

 * fasta format
 * separator format : sequence name (the separator can be a space or a tabulation or a comma or a semi-comma)
 * raw format : only the tag/kmer sequence

You must use option '-i' to specify the tag/kmer file.  You can provide the tag/kmer
file via the standard input by using '-i -' as filename.  If you file is
gziped, you can pass directly with the '-i mytags.gz' option or the pipe if
needed, but if it is in other compression format, uncompress the file with the
right tool and pass to countTags via the pipe and option '-i -'.

All tags/kmers must have at least the K-mer length, if too short, tags/kmers are discarded.
They are print to STDERR.

The maximum authorize tag length is 32bp (one integer).

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

