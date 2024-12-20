Changelog
=========

countTags/0.6.1 (2024-11-30) [Anthony Boureux]
----------------------------

## New

 * Output all reads matching one or more kmer in tabular and fasq format.

### Fix

 * Correct default k value to 31.

countTags/0.6.0 (2024-11-30) [Anthony Boureux]
----------------------------


### New

 * with option -r: output all data (header, sequence and quality code) to build a true fastq file
   - for now 

### Changes

 * Improve help and README.
 * Improve make test output.
  * better output with color
   * correct test with diff
 * Add code to debug, commented for now.

### Fix

 * Correct 20241113-01 but read smaller than k.
   - error with nline counter
   - corrected with bug 20241120-01
 * Correct bug 20241120-01 not checking truncated fastq.
 * Correct 20241113-01: forget to count the current line when read smaller than k.
 * Improve comments.
 * When output reads, output only once the read.
  * The read was output each time a kmer was found in the read
 * Minor: Update doc to install with zstr lib.

### Todo

 * Count twice when paired and overlapping.


countTags/0.5.2 (2020-04-19) [Anthony Boureux]
----------------------------


### New

 * Can use now ,; to separate tag and name.

### Changes

 * Improve color the make test.
 * Use namespace to avoid all std::
 * Minor: Correct g++ warnings.
 * Improve comments.

### Fix

 * 20200419-01: do not output number of read in summary.
   - closed: commit 459de402587.
 * Minor: Exit with code 10 when error reading tags or fastq.
 * Minor: Correct error with submodule.
 * Zstr: Integrated correctly the submodule.

countTags/0.5.1 (2020-03-26) [Anthony Boureux]
----------------------------


### Changes

 * Typo correction.
 * Remove require for --merge-counts when using --merge-counts-colname.
   - setting the second, set the first

### Fix

 * 20200326-01 : No normalization (closed).

countTags/0.5.0 (2020-03-10) [Anthony Boureux]
----------------------------


### Changes

 * To get statistic, use -v or --summary, remove other output to stdin other than table count.
 * Remove dirname if present in fastq filename.


countTags/0.4.7 (2020-03-10) [Anthony Boureux]
----------------------------


### New

 * Manage Tags in gz or not format, also from STDIN.
 * Add option to analyze all kmers in tags provided.
 * Use option -o to output merged tsv files in a file, otherwise output to stdout.
 * Add option to group summary like tsv files.
 * Shell script to merge count of countTags in one file.

### Fix

 * Improve help for new gzip format for TAGS.
 * Correct g++ warnings.
 * Do not need -s option, the script take care of the filename.
 * Add option nostranded, usefull in script.
  * this option exist only to avoid an empty option when testing multiple options

### Other

 * Add submodule info.
 * Minor: Improve help and verbose mode.
 * Minor: fix error in usage help.
 * Add return code for function.


countTags/0.4.5 (2020-03-08) [Anthony Boureux]
----------------------------

 * Add option BILLION to normalyze.
 * Remove short option '-m', keep only long format.
 * Update the documentation.
 * Add normalization by KBP.
 * Improve help messages.
 * [New] Output version number in summary.

countTags/0.4.3 (2019-11-22) [Anthony Boureux]
----------------------------

 * Add version number to help.

countTags/0.4.2 (2019-11-22) [Anthony Boureux]
----------------------------

 * New option summary
 * Improve code, help.
 * [Bug] mergeTagCounts do not manage tag name in input file.

countTags/0.4 (2019-11-12) [Anthony Boureux]
--------------------------

 * [New] Manage paired mode rf fr and ff.
   - work only with 2 fastq.
   * [Update] fastq files as to be name _1.fastq and _2.fastq in paired mode.
     only way to be sure to find _1.fastq and not ..._1...
 * Improve help message.
 * [Update] Add long tag names to check error.
 * [Change] Improve error output to stderr instead of stdout.
 * [New] Add tag with long header text.
 * [Change] Add help for us to write arguments list.
 * [Change] Output information only for debug.
 * [Change] Output total factor in stderr, not in the file.

countTags/0.3.5 (2019-08-27) [Anthony Boureux]
----------------------------

 * Minor correction, version 0.3.5.
 * Use more explicit name for filename to output reads.
 * Print column name when merging only if verbose.
 * Add Tags list for quality.
 * Correct date for AB copyright.
 * Add GPL copyrigth to each files.
 * New column merge format.
 * [Bug] Correction for help page.
 * [New] Add option to output a column name when merge-counts option is used.
 * [Bug] Resolved bug 20180508-02.

countTags/0.3.4 (2018-07-02) [Anthony Boureux]
----------------------------

 * [Bug] 20180625-01 : do not manage paired fastq file in stranded mode.
   - the first pair is usually in Reverse strand : Reverse-Forward
     But in this case countTags do not reverse the first pair, so it
     can not find a tag in this file.
 * [Bug] resolved 20180508-02: use stdin only if '-'
 * [Bug] resolved Bug 20180508-01.
 * [Bug] Do not check if no tag given.

countTags/0.3.3 (2018-05-07) [Anthony Boureux]
----------------------------

 * Add test with stranded option.
 * Generate temporary file in /tmp during test.
 * Do all tests in makefile.
 * Improve doc.
 * [Bug] test presence fastq file (#20180328-01)
 * [Bug] Should use pclose instead of fclose.
   - The fastq file is open via popen, should close it with pclose
   - Appenrently, do not generate an error or segfault
 * [New] Add option to store --reads matching kmer.
 * [Bug] 20180329-01, 20180326-01: Resolved, can use -h or -V without error.
 * [Bugs] conflict between -i and -h|-V.

countTags/0.3.2 (2018-03-29) [Anthony Boureux]
----------------------------

 * Can take from stdin
 * [Bug] Output on STDERR tag shorter than Kmer.
 * Correct the code in option::parser to make it working
 * Add a file to list Bugs.
 * Add option -i to check if a tag file is present.
 * On error, should return <> 0.
 * [Bug] We do not support kmer > 32.
 * v0.2 does not manage tag shorter than kmer
 * [Bug] Add too much sequences in fasta format.
 * First Test files to improve version checking.
 * [Bug] Miss a KMER_LENGTH in last merge.
 * Keep only tag with size > tag_length.
 * Improve help in README.
 * Sum count for all fastq files.
 * Join multiple countTags results in one file.

countTags/0.3 (2018-03-27) [Anthony Boureux]
--------------------------

 * Reduce if test.
 * Update optionparser.h from upstream to vers. 1.7.
 * Really improve program usage.
 * Can read tags from STDIN or a file.
 * Improve usage doc.
 * [Bug] fasta name not set properly.
 * [Bug] Correct the tag number according to the line number.
 * Add a default name for tag if require by print tag-names option.
 * Add verbose option and debug print.
 * Read fasta, raw or tsv file format for tags.
 * Stop on option error.
 * Add -v|--version option.
 * Add -h option for help.
 * Say error on bad line arguments.
 * Remove unused options and update documentation.
 * First standalone version of countTags.

