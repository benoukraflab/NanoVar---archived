AUTOMAKE_OPTIONS = foreign

DIR := ${CURDIR}

hsblastdir = $(DIR)/blast/hs-blastn-src

if LINUX
xmd5sum = "${md5sumcheck} -t"
endif
if DARWIN
xmd5sum = "${md5sumcheck} -r"
endif
.PHONY: buildblast pythonvirten scriptupdate movin clean distclean $(DIR) $(hsblastdir) $(SUBDIRS) check checkversions
SUBDIRS = scripts

all: nnvr
	@cp $(DIR)/scripts/nanovar $(DIR)/
	@echo "Compilation completed"

nnvr: checkversions buildblast pythonvirten movin scriptupdate

checkversions:
	@echo ${samversion} | awk -F'.' '{if ($$1<1 || ($$1==1 && $$2<4)) {print "Dependency error: samtools required to be version 1.4 or greater"; exit 1}}'
	@echo ${bedversion} | awk -F'.' '{if ($$1<2 || ($$1==2 && $$2<26)) {print "Dependency error: bedtools required to be version 2.26.0 or greater"; exit 1}}'
if DISABLE_BOWTIE2
	@echo "Compiling with mode --disable-bowtie2"
else
	@echo ${bowversion} | awk -F'.' '{if ($$1<2 || ($$1==2 && $$2<3) || ($$1==2 && $$2==3 && $$3<4)) {print "Dependency error: bowtie2 required to be version 2.3.4 or greater"; exit 1}}'
endif

buildblast: $(hsblastdir)
	@echo "Building blast..."
	@$(MAKE) -C $(hsblastdir)
	@echo "blast finished!"

pythonvirten:
	@echo "Generating python 2.7 virtual environment"
	$(DIR)/virtualenv-15.2.0/virtualenv.py --python=python2.7 $(DIR)/nv_virtualenv
	@source $(DIR)/nv_virtualenv/bin/activate; \
	pip install --upgrade pip; \
	pip install h5py; \
	pip install backports.weakref==1.0rc1; \
	pip install --upgrade tensorflow; \
	pip install keras; \
	pip install matplotlib; \
	pip install seaborn; \
	pip install statistics; \
	pip install scipy; \
	deactivate
	@echo "Finished installing virtual environment"

movin:
if LINUX
	@zcat $(DIR)/blast_bin/makeblastdb-linux.gz > $(DIR)/blast/makeblastdb
	@zcat $(DIR)/blast_bin/windowmasker-linux.gz > $(DIR)/blast/windowmasker
	@echo "Moved-linux"
endif
if DARWIN
	@zcat $(DIR)/blast_bin/makeblastdb-mac.gz > $(DIR)/blast/makeblastdb
	@zcat $(DIR)/blast_bin/windowmasker-mac.gz > $(DIR)/blast/windowmasker
	@echo "Moved-mac"	
endif
	@chmod +x $(DIR)/blast/makeblastdb
	@chmod +x $(DIR)/blast/windowmasker

scriptupdate: ${SUBDIRS}
	@echo "Updating NanoVar..."
	@$(MAKE) -C $(SUBDIRS)
	@echo "NanoVar updated"

check:
	@echo "Testing unit..."
	@zcat $(DIR)/test/ref_test.fa.gz > $(DIR)/test/ref.fa
	@echo "Pre-testing..."
	@cp -r $(DIR)/test/pre-check-files $(DIR)/test/pre-check
	@${samtoolcheck} faidx $(DIR)/test/ref.fa || (echo "Error: samtools faidx test failed, please ensure samtools is functional" && exit 1)
	@$(DIR)/blast/makeblastdb -in $(DIR)/test/pre-check/correct.fa -input_type fasta -dbtype nucl -logfile $(DIR)/test/pre-check/tmp || (echo "Error: makeblastdb binaries test failed, please download Blast v2.7.1+ source, build it and copy 'makeblastdb' binaries into $(DIR)/blast/" && exit 1)
	@$(DIR)/blast/windowmasker -in $(DIR)/test/pre-check/correct.fa -infmt blastdb -mk_counts -out $(DIR)/test/pre-check/correct.counts 2>$(DIR)/test/pre-check/tmp || (echo "Error: windowmasker binaries test1 failed, please download Blast v2.7.1+ source, build it and copy 'windowmasker' binaries into $(DIR)/blast/" && exit 1)
	@$(DIR)/blast/windowmasker -in $(DIR)/test/pre-check/correct.counts -sformat obinary -out $(DIR)/test/pre-check/correct.counts.obinary -convert 2>$(DIR)/test/pre-check/tmp 2>$(DIR)/test/pre-check/tmp || (echo "Error: windowmasker binaries test2 failed, please download Blast v2.7.1+ source, build it and copy 'windowmasker' binaries into $(DIR)/blast/" && exit 1)
if DISABLE_BOWTIE2
	@cat $(DIR)/test/pre-check/sam | ${samtoolcheck} view -@ 1 -Sb -f 2 - -o $(DIR)/test/pre-check/bam || (echo "Error: samtools view test failed, please ensure samtools is functional" && exit 1)
	@cat $(DIR)/test/pre-check/bam | ${samtoolcheck} sort -@ 1 -n - -o $(DIR)/test/pre-check/bam.sort || (echo "Error: samtools sort test failed, please ensure samtools is functional" && exit 1)
	@${samtoolcheck} fastq -@ 1 -1 $(DIR)/test/pre-check/bam.sort.fq1 -2 $(DIR)/test/pre-check/bam.sort.fq2 $(DIR)/test/pre-check/bam.sort 2>$(DIR)/test/pre-check/tmp || (echo "Error: samtools fastq test failed, please ensure samtools is functional" && exit 1)
	@${bedtoolcheck} genomecov -bg -pc -ibam $(DIR)/test/pre-check/bam > $(DIR)/test/pre-check/bedgraph || (echo "Error: bedtools genomecov test failed, please ensure bedtools is functional" && exit 1)
	@cat $(DIR)/test/pre-check/bedgraph | ${bedtoolcheck} sort -i - > $(DIR)/test/pre-check/bedgraph.sort || (echo "Error: bedtools sort test failed, please ensure bedtools is functional" && exit 1)
	@${bedtoolcheck} map -a $(DIR)/test/pre-check/overlap.ANN0.6.tsv -b $(DIR)/test/pre-check/bedgraph.sort -c 4 -o mean -null 0 > $(DIR)/test/pre-check/cov.tsv || (echo "Error: bedtools map test failed, please ensure bedtools is functional" && exit 1)
	@cat $(DIR)/test/pre-check/cov.tsv | awk -F'\t' "{if (0<=$$14 && $$14<=100) print $$0}" > $(DIR)/test/pre-check/cov.tsv2 || (echo "Error: awk test1 failed, please ensure awk is functional" && exit 1)
	@cat $(DIR)/test/pre-check/cov.tsv2 | sort -t $$'\t' -k9,9 -V > $(DIR)/test/pre-check/cov.tsv3 || (echo "Error: sort test failed, please ensure 'sort' is functional" && exit 1)
	@${bedtoolcheck} intersect -wa -wb -a $(DIR)/test/pre-check/bed1.sort -b $(DIR)/test/pre-check/bed2.sort > $(DIR)/test/pre-check/intersect || (echo "Error: bedtools intersect test failed, please ensure bedtools is functional" && exit 1)
	@cat $(DIR)/test/pre-check/intersect | awk -F'\t' '{if ($$4 != $$9) print $$0}' | awk -F'\t' '{if ($$5 == $$10) print $$0}' > $(DIR)/test/pre-check/intersect2 || (echo "Error: awk test2 failed, please ensure awk is functional" && exit 1)
	@cat $(DIR)/test/pre-check/cov.tsv3 | awk -F'\t' "{print $$0\"\t.\"}" > $(DIR)/test/pre-check/cov.tsv4 || (echo "Error: awk test3 failed, please ensure awk is functional" && exit 1)
	@trial=0.05; \
	echo $$trial | awk '{print 1 - $$1}' > $(DIR)/test/pre-check/trial || (echo "Error: awk test4 failed, please ensure awk is functional" && exit 1)
	@rm -r $(DIR)/test/pre-check
	@echo "Pre-testing passed"
	@echo "Post-testing..."
	@$(DIR)/scripts/nanovar -r $(DIR)/test/ref.fa -l $(DIR)/test/lread_test.fa.gz -o $(DIR)/test/testout -w 200 || (echo "Error: Intermediate test files corrupted, try 'make distclean' and rebuild" && exit 1)
	@$xmd5sum $(DIR)/test/ref.fa $(DIR)/test/lread_test.fa.gz $(DIR)/test/sread1.names_48k.fa.gz $(DIR)/test/sread2.names_48k.fa.gz $(DIR)/test/ref.fa.fai $(DIR)/test/ref.fa.nhr $(DIR)/test/ref.fa.nin $(DIR)/test/ref.fa.nsq $(DIR)/test/ref.fa.bwt $(DIR)/test/ref.fa.header $(DIR)/test/ref.fa.sa $(DIR)/test/ref.fa.sequence $(DIR)/test/testout/nanovar_run/hsblast_longreads/ref.fa.counts $(DIR)/test/testout/nanovar_run/hsblast_longreads/ref.fa.counts.obinary $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref.SV0.95.parse150.tsv $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref.SV0.95.parse150.1-13.7.overlap.tsv $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref.SV0.95.parse150.1-13.7.overlap.ANN.tsv $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref.SV0.95.parse150.1-13.7.overlap.ANN0.tsv $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref.SV0.95.parse150.1-13.7.overlap.ANN0.cov.tsv > $(DIR)/test/md5.out
	@tail -n +9 $(DIR)/test/testout/nanovar_results/lread_test.output.total.vcf | $xmd5sum >> $(DIR)/test/md5.out
	@tail -n +9 $(DIR)/test/testout/nanovar_results/lread_test.output.filtered-2.6.vcf | $xmd5sum >> $(DIR)/test/md5.out
	@$xmd5sum $(DIR)/test/testout/nanovar_results/lread_test.svread-overlap.tsv >> $(DIR)/test/md5.out
	@tail -n +64 $(DIR)/test/testout/nanovar_results/lread_test.output.report.html | $xmd5sum >> $(DIR)/test/md5.out
	@$xmd5sum $(DIR)/test/testout/nanovar_results/figures/depth_of_coverage.png $(DIR)/test/testout/nanovar_results/figures/read_length_dist.png $(DIR)/test/testout/nanovar_results/figures/scatter1.png $(DIR)/test/testout/nanovar_results/figures/scatter2.png $(DIR)/test/testout/nanovar_results/figures/sv_lengths.png $(DIR)/test/testout/nanovar_results/figures/sv_type_donut.png >> $(DIR)/test/md5.out
if LINUX
	@${pythoncheck} $(DIR)/scripts/nv_test.py $(DIR)/test/md5.out b0 $(DIR) linux
endif
if DARWIN
	@${pythoncheck} $(DIR)/scripts/nv_test.py $(DIR)/test/md5.out b0 $(DIR) mac
endif
	@if [ -s $(DIR)/test/AllPassed ]; then \
	echo "All passed!"; \
	else \
	echo "Error: Make check failed, try 'make distclean' and re-build NanoVar" && exit 1; \
	fi
	@rm $(DIR)/test/AllPassed
	@rm -r $(DIR)/test/testout
	@rm $(DIR)/test/ref.fa*
else
	@${bowtie2buildcheck} --threads 1 -q $(DIR)/test/pre-check/correct.fa $(DIR)/test/pre-check/correct || (echo "Error: bowtie2-build test failed, please ensure bowtie2-build is functional" && exit 1)
	@${bowtie2check} -p 1 -q --no-unal --quiet -x $(DIR)/test/pre-check/correct -1 $(DIR)/test/pre-check/fastq1 -2 $(DIR)/test/pre-check/fastq2 > $(DIR)/test/pre-check/sam || (echo "Error: bowtie2 test failed, please ensure bowtie2 is functional" && exit 1)
	@cat $(DIR)/test/pre-check/sam | ${samtoolcheck} view -@ 1 -Sb -f 2 - -o $(DIR)/test/pre-check/bam || (echo "Error: samtools view test failed, please ensure samtools is functional" && exit 1)
	@cat $(DIR)/test/pre-check/bam | ${samtoolcheck} sort -@ 1 -n - -o $(DIR)/test/pre-check/bam.sort || (echo "Error: samtools sort test failed, please ensure samtools is functional" && exit 1)
	@${samtoolcheck} fastq -@ 1 -1 $(DIR)/test/pre-check/bam.sort.fq1 -2 $(DIR)/test/pre-check/bam.sort.fq2 $(DIR)/test/pre-check/bam.sort 2>$(DIR)/test/pre-check/tmp || (echo "Error: samtools fastq test failed, please ensure samtools is functional" && exit 1)
	@${bedtoolcheck} genomecov -bg -pc -ibam $(DIR)/test/pre-check/bam > $(DIR)/test/pre-check/bedgraph || (echo "Error: bedtools genomecov test failed, please ensure bedtools is functional" && exit 1)
	@cat $(DIR)/test/pre-check/bedgraph | ${bedtoolcheck} sort -i - > $(DIR)/test/pre-check/bedgraph.sort || (echo "Error: bedtools sort test failed, please ensure bedtools is functional" && exit 1)
	@${bedtoolcheck} map -a $(DIR)/test/pre-check/overlap.ANN0.6.tsv -b $(DIR)/test/pre-check/bedgraph.sort -c 4 -o mean -null 0 > $(DIR)/test/pre-check/cov.tsv || (echo "Error: bedtools map test failed, please ensure bedtools is functional" && exit 1)
	@cat $(DIR)/test/pre-check/cov.tsv | awk -F'\t' "{if (0<=$$14 && $$14<=100) print $$0}" > $(DIR)/test/pre-check/cov.tsv2 || (echo "Error: awk test1 failed, please ensure awk is functional" && exit 1)
	@cat $(DIR)/test/pre-check/cov.tsv2 | sort -t $$'\t' -k9,9 -V > $(DIR)/test/pre-check/cov.tsv3 || (echo "Error: sort test failed, please ensure 'sort' is functional" && exit 1)
	@${bedtoolcheck} intersect -wa -wb -a $(DIR)/test/pre-check/bed1.sort -b $(DIR)/test/pre-check/bed2.sort > $(DIR)/test/pre-check/intersect || (echo "Error: bedtools intersect test failed, please ensure bedtools is functional" && exit 1)
	@cat $(DIR)/test/pre-check/intersect | awk -F'\t' '{if ($$4 != $$9) print $$0}' | awk -F'\t' '{if ($$5 == $$10) print $$0}' > $(DIR)/test/pre-check/intersect2 || (echo "Error: awk test2 failed, please ensure awk is functional" && exit 1)
	@cat $(DIR)/test/pre-check/cov.tsv3 | awk -F'\t' "{print $$0\"\t.\"}" > $(DIR)/test/pre-check/cov.tsv4 || (echo "Error: awk test3 failed, please ensure awk is functional" && exit 1)
	@trial=0.05; \
	echo $$trial | awk '{print 1 - $$1}' > $(DIR)/test/pre-check/trial || (echo "Error: awk test4 failed, please ensure awk is functional" && exit 1)
	@rm -r $(DIR)/test/pre-check
	@echo "Pre-testing passed"
	@echo "Post-testing..."
	@$(DIR)/scripts/nanovar -r $(DIR)/test/ref.fa -l $(DIR)/test/lread_test.fa.gz -s1 $(DIR)/test/sread1.names_48k.fa.gz -s2 $(DIR)/test/sread2.names_48k.fa.gz -o $(DIR)/test/testout -w 200 || (echo "Error: Intermediate test files corrupted, try 'make distclean' and rebuild" && exit 1)
	@$xmd5sum $(DIR)/test/ref.fa $(DIR)/test/lread_test.fa.gz $(DIR)/test/sread1.names_48k.fa.gz $(DIR)/test/sread2.names_48k.fa.gz $(DIR)/test/ref.fa.fai $(DIR)/test/ref.fa.nhr $(DIR)/test/ref.fa.nin $(DIR)/test/ref.fa.nsq $(DIR)/test/ref.fa.bwt $(DIR)/test/ref.fa.header $(DIR)/test/ref.fa.sa $(DIR)/test/ref.fa.sequence $(DIR)/test/testout/nanovar_run/hsblast_longreads/ref.fa.counts $(DIR)/test/testout/nanovar_run/hsblast_longreads/ref.fa.counts.obinary $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref.SV0.95.parse150.tsv $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref.SV0.95.parse150.1-13.7.overlap.tsv $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref.SV0.95.parse150.1-13.7.overlap.ANN.tsv $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref.SV0.95.parse150.1-13.7.overlap.ANN0.tsv $(DIR)/test/testout/nanovar_run/bowtie2_shortreads/lread_test.overlap.ANN0.correct.1.bt2 $(DIR)/test/testout/nanovar_run/bowtie2_shortreads/lread_test.overlap.ANN0.correct.2.bt2 $(DIR)/test/testout/nanovar_run/bowtie2_shortreads/lread_test.overlap.ANN0.correct.3.bt2 $(DIR)/test/testout/nanovar_run/bowtie2_shortreads/lread_test.overlap.ANN0.correct.4.bt2 $(DIR)/test/testout/nanovar_run/bowtie2_shortreads/lread_test.overlap.ANN0.correct.rev.1.bt2 $(DIR)/test/testout/nanovar_run/bowtie2_shortreads/lread_test.overlap.ANN0.correct.rev.2.bt2 > $(DIR)/test/md5.out
	@${samtoolcheck} view $(DIR)/test/testout/nanovar_run/bowtie2_shortreads/lread_test.overlap.ANN0.correct.bow-longr-f2-ns.bam | $xmd5sum >>$(DIR)/test/md5.out
	@$xmd5sum $(DIR)/test/testout/nanovar_run/bowtie2_shortreads/lread_test.overlap.ANN0.correct.bow-longr-f2-ns_1.fastq $(DIR)/test/testout/nanovar_run/bowtie2_shortreads/lread_test.overlap.ANN0.correct.bow-longr-f2-ns_2.fastq $(DIR)/test/testout/nanovar_run/hsblast_longreads/lread_test.hsblast-ref.SV0.95.parse150.1-13.7.overlap.ANN0.0-100.tsv >> $(DIR)/test/md5.out
	@tail -n +9 $(DIR)/test/testout/nanovar_results/lread_test.output.total.vcf | $xmd5sum >> $(DIR)/test/md5.out
	@tail -n +9 $(DIR)/test/testout/nanovar_results/lread_test.output.filtered-2.6.vcf | $xmd5sum >> $(DIR)/test/md5.out
	@$xmd5sum $(DIR)/test/testout/nanovar_results/lread_test.svread-overlap.tsv >> $(DIR)/test/md5.out
	@tail -n +64 $(DIR)/test/testout/nanovar_results/lread_test.output.report.html | $xmd5sum >> $(DIR)/test/md5.out
	@$xmd5sum $(DIR)/test/testout/nanovar_results/figures/depth_of_coverage.png $(DIR)/test/testout/nanovar_results/figures/read_length_dist.png $(DIR)/test/testout/nanovar_results/figures/scatter1.png $(DIR)/test/testout/nanovar_results/figures/scatter2.png $(DIR)/test/testout/nanovar_results/figures/sv_lengths.png $(DIR)/test/testout/nanovar_results/figures/sv_type_donut.png >> $(DIR)/test/md5.out
if LINUX
	@${pythoncheck} $(DIR)/scripts/nv_test.py $(DIR)/test/md5.out b1 $(DIR) linux
endif
if DARWIN
	@${pythoncheck} $(DIR)/scripts/nv_test.py $(DIR)/test/md5.out b1 $(DIR) mac
endif
	@if [ -s $(DIR)/test/AllPassed ]; then \
	echo "All passed!"; \
	else \
	echo "Error: Make check failed, try 'make distclean' and re-build NanoVar" && exit 1; \
	fi
	@rm $(DIR)/test/AllPassed
	@rm -r $(DIR)/test/testout
	@rm $(DIR)/test/ref.fa*
endif
 
clean: $(SUBDIRS) $(hsblastdir)
	@$(MAKE) -C $(SUBDIRS) clean
	@$(MAKE) -C $(hsblastdir) clean
	@rm -rf $(DIR)/nv_virtualenv
	@rm -f $(DIR)/nanovar
	@rm -rf $(DIR)/test/testout
	@rm -f $(DIR)/test/ref.fa*
	@rm -f $(DIR)/test/AllPassed
	@rm -rf $(DIR)/test/pre-check
	@echo "CLEAN"

distclean: $(SUBDIRS) $(hsblastdir)
	@$(MAKE) -C $(SUBDIRS) distclean
	@$(MAKE) -C $(hsblastdir) distclean
	@rm -rf $(DIR)/nv_virtualenv
	@rm -f $(DIR)/nanovar
	@rm -rf $(DIR)/test/testout
	@rm -f $(DIR)/test/ref.fa*
	@rm -f $(DIR)/test/AllPassed
	@rm -rf $(DIR)/test/pre-check
	@rm -f $(DIR)/blast/makeblastdb $(DIR)/blast/windowmasker 
	@rm -f $(DIR)/Makefile
	@rm -f $(DIR)/test/md5.out
	@rm -f $(DIR)/config.log
	@rm -rf $(DIR)/autom4te.cache
	@rm -f $(DIR)/config.status
	@echo "DISTCLEAN"
