bin_SCRIPTS = nanovar
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = nanovar.in

MDIR = ${abs_top_builddir}

do_sub_wout_bow = @sed -e 's|bowtie2path=LOADING|bowtie2path=$(MDIR)/scripts/bowtie2_dummy.py|g' \
	-e 's|bowtie2build=LOADING|bowtie2build=$(MDIR)/scripts/bowtie2_dummy.py|g' \
	-e 's|bow=0|bow=1|g' \
	-e 's|Optional|Disabled|g' \
	-e 's|maindir=LOADING|maindir=$(MDIR)|g' \
	-e 's|samtoolpath=LOADING|samtoolpath=${samtoolcheck}|g' \
	-e 's|bedtoolpath=LOADING|bedtoolpath=${bedtoolcheck}|g'

do_sub_with_bow = @sed -e 's|bowtie2path=LOADING|bowtie2path=$(bowtie2check)|g' \
        -e 's|bowtie2build=LOADING|bowtie2build=$(bowtie2buildcheck)|g' \
        -e 's|maindir=LOADING|maindir=$(MDIR)|g' \
        -e 's|samtoolpath=LOADING|samtoolpath=${samtoolcheck}|g' \
        -e 's|bedtoolpath=LOADING|bedtoolpath=${bedtoolcheck}|g'

nanovar: nanovar.in Makefile
if DISABLE_BOWTIE2
	$(do_sub_wout_bow) < ${CURDIR}/nanovar.in > nanovar
	@chmod +x nanovar
else
	$(do_sub_with_bow) < ${CURDIR}/nanovar.in > nanovar
	@chmod +x nanovar
endif
