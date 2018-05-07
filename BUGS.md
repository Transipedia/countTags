# [X] 20180329-02 : No information is output when tag are drop because they are shorter than kmer
 * Output on sdterr, tag shorter than kmer (commit 67546ef9d1ccf).

# [X] 20180329-01 : When asking for -h -V, got error for nothing with -i
 * due to Bug #20180326-01 (commit 451314ffc4a).

# [ ] 20180328-01 : Do not Test if fastq.gz files are present

# [ ] 20180328-02 : Get the shorter alphabetic tag in the output
 * Solution 1:
    - use 1 bit to store if the tag is forward or reverse-complement
    - so only 31 bits for the tag.

# [-] 20180326-02 : When kmer length > 32, sequences are scramble
 * at this time add a message and exit if kmer > 32 (commit bf8e7f80105)
 * Problem probably with storrage UINT32, as pointed by JA in Readme

# [X] 20180326-01 : Do not test if there is a tags file, always take first argument as tag file
 * Solution: put the tag filename as argument with option '-i'
 * Not Working as expected: the Arg::Required is not working
 * So for now, test with a 'if' condition, in main code (commit 451314ffc4a).
