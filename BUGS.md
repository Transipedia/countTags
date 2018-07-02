# [ ] 20180625-01 : Do not manage paired fastq files in stranded mode
 * Paired fastq file are not considered, always using forward strand for the two pair

# [X] 20180508-01 : Do not stop if no tag given
 * Add a test to exit if no tag given in tag file or stdin (commit 5dd253e)

# [X] 20180508-02 : Try to open stdin if tag filename is oneletter
 * Test filename is '-' if stdin, and not filename size =1 (commit 036e0192)

# [X] 20180329-02 : No information is output when tag are drop because they are shorter than kmer
 * Output on sdterr, tag shorter than kmer (commit 67546ef9d1ccf).

# [X] 20180329-01 : When asking for -h -V, got error for nothing with -i
 * due to Bug #20180326-01 (commit 451314ffc4a).

# [X] 20180328-01 : Do not Test if fastq.gz files are present
 * done in commit 63b14cd3eb.

# [ ] 20180328-02 : Get the shorter alphabetic tag in the output
 * Solution 1:
    - use 1 bit to store if the tag is forward or reverse-complement
    - but only 31 bits for the tag.
 * Solution 2:
    - store in hash: DNAid -> sequence
    - store in hash, or array if get tag line number: 0 for F, 1 for R, and reverse or not the DNAid if required

# [-] 20180326-02 : When kmer length > 32, sequences are scramble
 * at this time add a message and exit if kmer > 32 (commit bf8e7f80105)
 * Problem probably with storrage UINT32, as pointed by JA in Readme

# [X] 20180326-01 : Do not test if there is a tags file, always take first argument as tag file
 * Solution: put the tag filename as argument with option '-i'
 * Not Working as expected: the Arg::Required is not working
 * So for now, test with a 'if' condition, in main code (commit 451314ffc4a).
