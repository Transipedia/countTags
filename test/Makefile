include Makefile.ref
SHELL:=/bin/bash

all: clean test_ok test_stranded test_error_parsing test_error_fastq test_reads
	@if [ -e /tmp/$(ofile).error ]; then echo -e "$(RED)ERROR HAPPENS$(NoColor)" ; else echo -e "\n$(GREEN)ALL FINE$(NoColor)" ; rm /tmp/$(ofile)*; fi

test_ok:
	@echo -e "$(PURPLE)Test everything is fine:$(NoColor)"
	@echo -e "$(PURPLE)# test with tag-names$(NoColor)"
	@for f in fas tsv csv; do ../bin/countTags --tag-names -k 30 -i $(okfile).$$f test_?.fastq.gz 2>/tmp/$(ofile)_names_$${f}_ok.cerr | sort > /tmp/$(ofile)_names_$${f}_ok.txt; done
	@for i in /tmp/$(ofile)_names_*_ok.txt; do echo -e "$(BROWNBG)$$i$(NoColor):"; diff $$i ok/$(ofile)_with_names_ok.txt; if [ $$? -eq 0 ]; then echo -e "$(GREEN)ok$(NoColor)"; else echo -e "$(RED)ERROR$(NoColor)"; echo "1" > /tmp/$(ofile).error; fi;done
	@echo -e "$(PURPLE)# test without tag-names$(NoColor)"
	@for f in fas tsv csv raw; do ../bin/countTags -k 30 -i $(okfile).$$f test_?.fastq.gz 2>/tmp/$(ofile)_without_$${f}_ok.cerr | sort > /tmp/$(ofile)_without_$${f}_ok.txt; done
	@for i in /tmp/$(ofile)_without_*.txt; do echo -e "$(BROWNBG)$$i$(NoColor):"; diff $$i ok/$(ofile)_without_names_ok.txt; if [ $$? -eq 0 ]; then echo -e "$(GREEN)ok$(NoColor)"; else echo -e "$(RED)ERROR$(NoColor)"; echo "1" > /tmp/$(ofile).error; fi;done

test_error_parsing:
	@echo -e "$(PURPLE)Test error in parsing argument:$(NoColor)"
	@../bin/countTags --tags_names --tag-format 2 -k 30 -m 10 $(okfile).raw test_?.fastq.gz 2>/dev/null >/dev/null ; if [ $$? -eq 1 ]; then echo -e "$(GREEN)ok$(NoColor)"; else echo -e "$(RED)ERROR$(NoColor)"; fi

test_error_fastq:
	@echo -e "$(PURPLE)Test error with trunc fastq:$(NoColor)"
	@for i in 2 3 4; do echo -e "$(BROWNBG)Missing line $$i in fastq$(NoColor):"; ../bin/countTags -k 30 -i $(okfile).raw test_trunc$$i.fastq.gz 2>/dev/null >/dev/null ; if [ $$? -eq $$i ]; then echo -e "$(GREEN)ok$(NoColor)"; else echo -e "$(RED)ERROR$(NoColor)"; echo "1" > /tmp/$(ofile).error; fi; done

test_stranded:
	@echo -e "$(PURPLE)Test with stranded option:$(NoColor)"
	@for f in fas tsv csv; do ../bin/countTags --tag-names -k 30 --stranded -i $(okfile).$$f test_?.fastq.gz 2>/tmp/$(ofile)_names_$${f}_okstranded.cerr | sort > /tmp/$(ofile)_names_$${f}_okstranded.txt; done
	@for i in /tmp/$(ofile)_names_*_okstranded.txt; do echo -e "$(BROWNBG)$$i$(NoColor)"; diff $$i stranded/$(ofile)_with_names_ok.txt; if [ $$? -eq 0 ]; then echo -e "$(GREEN)ok$(NoColor)"; else echo -e "$(RED)ERROR$(NoColor)"; echo "1" > /tmp/$(ofile).error; fi;done

test_reads:
	@echo -e "$(PURPLE)Test reads output option:$(NoColor)"
	@../bin/countTags --tag-names -k 30 -i $(okfile).csv -r /tmp/$(ofile)-reads test_?.fastq.gz 2>/tmp/$(ofile)_names_$${f}_reads.cerr > /dev/null
	@for i in /tmp/$(ofile)-reads* ; do echo -e "$(BROWNBG)$$i$(NoColor):"; diff $$i ok/$$(basename $$i); if [ $$? -eq 0 ]; then echo -e "$(GREEN)ok$(NoColor)"; else echo -e "$(RED)ERROR$(NoColor)"; echo "1" > /tmp/$(ofile).error; fi;done

init:
	rm -f $(okfile).*
	jellyfish count -s 3M -m 30 -o $(jf) <(zcat test_*.fastq.gz)
	jellyfish stats $(jf)
	jellyfish dump -L 20 -c -o $(jfcount) $(jf)
	# Generate reference line $(okfile).csv
	# Add uniq tag number at the name end
	for i in $$(seq 1 15); do l=$$(echo $$(( ( RANDOM % $$(wc -l $(jfcount)| cut -d " " -f1) )  + 1 )) ); tail -n +$$l $(jfcount) | head -1; done | perl -p -e 'BEGIN {$$i=1};$$j=$$i++;s/\n/.$$j\n/' > $(okfile).csv
	# Add scrable read that will not match
	head -5 $(okfile).csv | tr '[ATGC]' '[GACT]' | sed 's/ / not/' >> $(okfile).csv
	# Add kmer with potential problems in header (created manually from $(okfile).csv
	cat TAGS_test_supl.csv >> $(okfile).csv

	# Add shorter tags than 30; should not be analyse
	tail -n +5 $(okfile).csv | head -4 | cut -c 1-15 | perl -p -e 'BEGIN {$$i=1};$$j=$$i++;s/\n/ shok.$$j\n/' >> $(okfile).csv
	# Add scramble shorter tags than 30; should not be analyse
	tail -n +5 $(okfile).csv | head -4 | cut -c 1-15 | tr 'A' 'G' | perl -p -e 'BEGIN {$$i=1};$$j=$$i++;s/\n/ shbad.$$j\n/' >> $(okfile).csv
	# Create tsv file format
	cat $(okfile).csv | sed 's/ /\t/' > $(okfile).tsv
	# Create fasta file format
	cat $(okfile).tsv | while read tag count; do echo -e ">$$count\n$$tag"; done > $(okfile).fas
	# Create raw file format
	cut -d ' ' -f1 $(okfile).csv > $(okfile).raw
	# Generate expected results
	../attic/$(countTagsReference) --tag-names -k 30 -i $(okfile).fas test_?.fastq.gz 2>ok/$(ofile)_with_names_ok.cerr | sort > ok/$(ofile)_with_names_ok.txt
	../attic/$(countTagsReference) -k 30 -i $(okfile).fas test_?.fastq.gz  2>ok/$(ofile)_without_names_ok.cerr | sort > ok/$(ofile)_without_names_ok.txt
	../attic/$(countTagsReference) --tag-names -k 30 -i $(okfile).fas --stranded test_?.fastq.gz  2>stranded/$(ofile)_with_names_ok.cerr | sort > stranded/$(ofile)_with_names_ok.txt

clean:
	@rm -f /tmp/$(ofile)*

.PHONY: clean
