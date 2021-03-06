include Makefile.ref

all: test_ok test_error_parsing test_stranded

test_ok:
	@echo "$(PURPLE)Test everything is fine:$(NoColor)"
	@echo "$(PURPLE)# test with tag-names$(NoColor)"
	@for f in fas tsv csv; do ../bin/countTags --tag-names -k 30 -i $(okfile).$$f test_?.fastq.gz | sort > /tmp/$(ofile)_names_$${f}_ok.txt; done
	for i in /tmp/$(ofile)_names_*_ok.txt; do echo "$(BROWNBG)$$i$(NoColor):"; diff $$i ok/$(ofile)_with_names_ok.txt; if [ $$? -eq 0 ]; then echo "$(GREEN)ok$(NoColor)"; else echo "$(RED)ERROR$(NoColor)"; fi;done
	@echo "$(PURPLE)# test without tag-names$(NoColor)"
	@for f in fas tsv csv raw; do ../bin/countTags -k 30 -i $(okfile).$$f test_?.fastq.gz | sort > /tmp/$(ofile)_without_$${f}_ok.txt; done
	for i in /tmp/$(ofile)_without_*.txt; do echo "$(BROWNBG)$$i$(NoColor):"; diff $$i ok/$(ofile)_without_names_ok.txt; if [ $$? -eq 0 ]; then echo "$(GREEN)ok$(NoColor)"; else echo "$(RED)ERROR$(NoColor)"; fi;done

test_error_parsing:
	@echo "$(PURPLE)Test error in parsing argument:$(NoColor)"
	../bin/countTags --tags_names --tag-format 2 -k 30 -m 10 $(okfile).raw test_?.fastq.gz

test_stranded:
	@echo "$(PURPLE)Test with stranded option:$(NoColor)"
	for f in fas tsv csv; do ../bin/countTags --tag-names -k 30 --stranded -i $(okfile).$$f test_?.fastq.gz | sort > /tmp/$(ofile)_names_$${f}_okstranded.txt; done
	for i in /tmp/$(ofile)_names_*_okstranded.txt; do echo "$(BROWNBG)$$i$(NoColor)"; diff $$i stranded/$(ofile)_with_names_ok.txt; if [ $$? -eq 0 ]; then echo "$(GREEN)ok$(NoColor)"; else echo "$(RED)ERROR$(NoColor)"; fi;done

init:
	rm -f $(okfile).*
	# use bash instead of sh to <()
	bash -c 'jellyfish count -s 3M -m 30 -o $(jf) <(zcat test_*.fastq.gz)'
	jellyfish stats $(jf)
	jellyfish dump -L 20 -c -o $(jfcount) $(jf)
	# Generate reference line $(okfile).csv
	# Add uniq tag number at the name end
	bash -c 'for i in $$(seq 1 15); do l=$$(echo $$(( ( RANDOM % $$(wc -l $(jfcount)| cut -d " " -f1) )  + 1 )) ); tail -n +$$l $(jfcount) | head -1; done' | perl -p -e 'BEGIN {$$i=1};$$j=$$i++;s/\n/.$$j\n/' > $(okfile).csv
	# Add scrable read that will not match
	head -5 $(okfile).csv | tr '[ATGC]' '[GACT]' | sed 's/ / not/' >> $(okfile).csv
	# Add kmer with potential problems in header (create manually from $(okfile).csv
	cat TAGS_test_supl.csv >> $(okfile).csv

	# Add shorter tags than 30; should not be analyse
	tail -n +5 $(okfile).csv | head -4 | cut -c 1-15 | perl -p -e 'BEGIN {$$i=1};$$j=$$i++;s/\n/ shok.$$j\n/' >> $(okfile).csv
	# Add scramble shorter tags than 30; should not be analyse
	tail -n +5 $(okfile).csv | head -4 | cut -c 1-15 | tr 'A' 'G' | perl -p -e 'BEGIN {$$i=1};$$j=$$i++;s/\n/ shbad.$$j\n/' >> $(okfile).csv
	# Create tsv file format
	cat $(okfile).csv | sed 's/ /\t/' > $(okfile).tsv
	# Create fasta file format
	cat $(okfile).tsv | while read tag count; do echo ">$$count\n$$tag"; done > $(okfile).fas
	# Create raw file format
	cut -d ' ' -f1 $(okfile).csv > $(okfile).raw
	# Generate expected results
	../attic/$(countTagsReference) --tag-names -k 30 -i $(okfile).fas test_?.fastq.gz | sort > ok/$(ofile)_with_names_ok.txt
	../attic/$(countTagsReference) -k 30 -i $(okfile).fas test_?.fastq.gz | sort > ok/$(ofile)_without_names_ok.txt
