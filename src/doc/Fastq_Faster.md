# Solution for the fastest way to read fastq gz files

## kseq
 * [website](http://attractivechaos.github.io/klib/#Kseq%3A%20stream%20buffer%20and%20FASTA%2FQ%20parser)
 * macro C

## fastq
 * fastq from stdin: [How To Efficiently Parse A Huge Fastq File](https://www.biostars.org/p/10353/)
 * [use a buffer](https://bytes.com/topic/c/answers/644744-fast-way-read-text-file-line-line)
 * (1) [use memory map buffer home made](http://create.stephan-brumme.com/portable-memory-mapping/) can be faster than wc, by mounting all file in memory
 * [Exist in boost also](https://stackoverflow.com/questions/10381610/using-boostiostreams-mapped-file-source-and-filtering-streambuf-to-decompress)

## multitreading
 * can use a omp pragram on a for loop to use multiple core [#1)
 * can use \<(zcat ) bash syntax to multithread the program

## Boost
 * [tuto](https://www.quantnet.com/threads/c-multithreading-in-boost.10028/)
 * [Boost Book](https://theboostcpplibraries.com/boost.thread-management)

## C++ pthread
 * [tuto](https://solarianprogrammer.com/2011/12/16/cpp-11-thread-tutorial/)

## grep
 * use LC_ALL=C to accelerate


