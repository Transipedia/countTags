include Makefile.ref

all: test_ok test_error_parsing

test_ok:
	# test with tag-names
	for f in fas tsv csv; do ../bin/countTags --tag-names -k 30 $(okfile).$$f test_?.fastq.gz > $(ofile)_names_$${f}_ok.txt; done
	for i in $(ofile)_names_*.txt; do echo $$i; diff $$i ok/$(ofile)_with_names_ok.txt;done
	# test without tag-names
	for f in fas tsv csv raw; do ../bin/countTags -k 30 $(okfile).$$f test_?.fastq.gz > $(ofile)_without_$${f}_ok.txt; done
	for i in $(ofile)_without_*.txt; do echo $i; diff $$i ok/$(ofile)_without_names_ok.txt;done

test_error_parsing:
	../bin/countTags --tags_names --tag-format 2 -k 30 -m 10 $(okfile).raw test_?.fastq.gz


init:
	rm -f $(okfile).*
	# use bash instead of sh to <()
	bash -c 'jellyfish count -s 3M -m 30 -o $(jf) <(zcat test_?.fastq.gz)'
	jellyfish stats $(jf)
	jellyfish dump -L 20 -c -o $(jfcount) $(jf)
	bash -c 'for i in $$(seq 1 15); do l=$$(echo $$(( ( RANDOM % $$(wc -l $(jfcount)| cut -d " " -f1) )  + 1 )) ); tail -n +$$l $(jfcount) | head -1; done' | perl -p -e 'BEGIN {$$i=1};$$j=$$i++;s/\n/.$$j\n/' > $(okfile).csv
	# Add uniq tag number at the name end
	# Add scrable read that will not match
	head -5 $(okfile).csv | tr '[ATGC]' '[GACT]' | sed 's/ / not/' >> $(okfile).csv

	# Generate good results, without shorter tag: HACK version 0.2 (TODO: to remove after reference is >=0.3)
	# Create fasta file format
	cat $(okfile).csv | sed 's/ /\t/' > $(okfile).tsv
	cat $(okfile).tsv | while read tag count; do echo ">$$count\n$$tag"; done > $(okfile).fas
	# Generate expected results
	../attic/$(countTagsReference) --tag-names -k 30 -m 100000 $(okfile).fas test_?.fastq.gz > ok/$(ofile)_with_names_ok.txt
	../attic/$(countTagsReference) -k 30 -m 100000 $(okfile).fas test_?.fastq.gz > ok/$(ofile)_without_names_ok.txt

	# Add shorter tags than 30; should not be analyse
	tail -n +5 $(okfile).csv | head -4 | cut -c 1-15 | perl -p -e 'BEGIN {$$i=1};$$j=$$i++;s/\n/ shok.$$j\n/'>> $(okfile).csv
	# Add scramble shorter tags than 30; should not be analyse
	tail -n +5 $(okfile).csv | head -4 | cut -c 1-15 | tr 'A' 'G' | perl -p -e 'BEGIN {$$i=1};$$j=$$i++;s/\n/ shbad.$$j\n/'>> $(okfile).csv
	# Create tsv file format
	cat $(okfile).csv | sed 's/ /\t/' > $(okfile).tsv
	# Create fasta file format
	cat $(okfile).tsv | while read tag count; do echo ">$$count\n$$tag"; done > $(okfile).fas
	# Create raw file format
	cut -d ' ' -f1 $(okfile).csv > $(okfile).raw
	# Generate expected results
	#../attic/$(countTagsReference) --tag-names -k 30 -m 100000 $(okfile).fas test_?.fastq.gz > ok/$(ofile)_with_names_ok.txt
	#../attic/$(countTagsReference) -k 30 -m 100000 $(okfile).fas test_?.fastq.gz > ok/$(ofile)_without_names_ok.txt
