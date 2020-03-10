# countTags - counting a set of k-mers into a set of FASTQ files

## Installation

1. Clone this repository : `git clone https://gitlab:get/counttags.git`
2. Compile the software : `make`
3. Place the binary somewhere acessible in your `$PATH`

## Usage

First you need to create a file containing the tags you want to quantify.
You can use three tags format:
 * fasta format
 * separator format : sequence name (the separator can be a space or a tabulation)
 * raw format : only the tag sequence

You must use option '-i' to specify the tag file.  You can provide the tags
file via the standard input by using '-i -' as filename.  If you file is
gziped, you can pass directly to the '-i mytags.gz' option or the pipe if
needed, but if it is in other compression format, uncompress the file with the
right tool and pass to countTags via the pipe and option '-i -'.

All tags must have at least the K-mer length, if too short, tags are discarded.
The maximum authorize tag length is 32bp.

K-mer length must be provided to countTags using the `-k INT` option.

For example :

`countTags -k 31 file.fa file1.fastq.gz file2.fastq.gz`

### Managing stranded files

By default, countTags count the canonical tag between the forward and the reverse tag
(the first one in alphabeticl order) and output this sequence in the result.

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

