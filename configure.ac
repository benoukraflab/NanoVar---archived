AC_INIT([NanoVar], [1.0.2], [chengyong.tham@u.nus.edu])
aclocal
automake --add-missing
AM_INIT_AUTOMAKE

#Test OS
AC_CANONICAL_HOST

os_linux=no
os_mac=no

case "${host_os}" in
	linux*)
		os_linux=yes
		;;
	darwin*)
		os_mac=yes
		;;
	*)
		AC_MSG_ERROR(["OS $host_os not supported"])
		;;
esac

AM_CONDITIONAL([LINUX], [test "x$os_linux" = xyes])
AM_CONDITIONAL([DARWIN], [test "x$os_mac" = xyes])

#gcc check
AC_PATH_PROG(gcccheck,gcc,no)
gccversion=$($gcccheck -dumpversion)
AS_IF([test -z "$gccversion"], [AC_MSG_ERROR([Error: Could not retrieve gcc version information])])
AC_SUBST([gccversion])

#ldd check
AC_PATH_PROG(lddcheck,ldd,no)
lddversion=$($lddcheck --version | head -n 1 | cut -d " " -f 3)
AS_IF([test -z "$lddversion"], [AC_MSG_ERROR([Error: Could not retrieve ldd version information])])
AC_SUBST([lddversion])

#Python2.7 check
AC_PATH_PROG(pythoncheck,python2.7,no)
AS_IF([test "x$pythoncheck" = xno], [AC_MSG_ERROR([Please ensure "python2.7" in \$PATH.])])

#seq check
AC_PATH_PROG(seqcheck,seq,no)
AS_IF([test "x$seqcheck" = xno], [AC_MSG_ERROR([Please ensure "seq" in \$PATH.])])

#touch check
AC_PATH_PROG(touchcheck,touch,no)
AS_IF([test "x$touchcheck" = xno], [AC_MSG_ERROR([Please ensure "touch" in \$PATH.])])
touch configure.ac aclocal.m4 configure Makefile.am Makefile.in  

#zcat check
AC_PATH_PROG(zcatcheck,zcat,no)
AS_IF([test "x$zcatcheck" = xno], [AC_MSG_ERROR([Please ensure "zcat" in \$PATH.])])

#md5sum check
AS_IF([test "x$os_linux" = xyes], [AC_PATH_PROG(md5sumcheck,md5sum,no)])
AS_IF([test "x$os_mac" = xyes], [AC_PATH_PROG(md5sumcheck,md5,no)])
AS_IF([test "x$md5sumcheck" = xno], [AC_MSG_ERROR([Please ensure "md5sum" or "md5" in \$PATH.])])

#samtools check
AC_ARG_WITH([samtools], AS_HELP_STRING([--with-samtools=/path/to/samtools], [Specify full path to samtools]), [echo Assigning samtools... "$with_samtools"], [AC_PATH_PROG(samtoolcheck,samtools,no)])
AS_IF([test ! -z "$with_samtools"], [AS_IF([test ! -x "$with_samtools" || test -d "$with_samtools"], [AC_MSG_ERROR([Error: samtools path does not point to an executable.])], [samtoolcheck=$with_samtools])])
AS_IF([test "x$samtoolcheck" = xno], [AC_MSG_ERROR([Error: samtools not found, please ensure "samtools" in \$PATH before installing NanoVar, or ./configure with --with-samtools=/path/to/samtools])])
samversion=$($samtoolcheck --version | head -n 1 | cut -d " " -f 2)
AS_IF([test -z "$samversion"], [AC_MSG_ERROR([Error: Could not retrieve samtools version information])])
AC_SUBST([samversion])

#bedtools check
AC_ARG_WITH([bedtools], AS_HELP_STRING([--with-bedtools=/path/to/bedtools], [Specify full path to bedtools]), [echo Assigning bedtools... "$with_bedtools"], [AC_PATH_PROG(bedtoolcheck,bedtools,no)])
AS_IF([test ! -z "$with_bedtools"], [AS_IF([test ! -x "$with_bedtools" || test -d "$with_bedtools"], [AC_MSG_ERROR([Error: bedtools path does not point to an executable.])], [bedtoolcheck=$with_bedtools])])
AS_IF([test "x$bedtoolcheck" = xno], [AC_MSG_ERROR([Error: bedtools not found, please ensure "bedtools" in \$PATH before installing NanoVar, or ./configure with --with-bedtools=/path/to/bedtools])])
bedversion=$($bedtoolcheck --version | sed 's|bedtools v||g')
AS_IF([test -z "$bedversion"], [AC_MSG_ERROR([Error: Could not retrieve bedtools version information])])
AC_SUBST([bedversion])

#bowtie2 check
AC_ARG_WITH([bowtie2], AS_HELP_STRING([--with-bowtie2=/path/to/bowtie2], [Specify full path to bowtie2]), [echo Assigning bowtie2... "$with_bowtie2"], [AC_PATH_PROG(bowtie2check,bowtie2,no)])
AS_IF([test ! -z "$with_bowtie2"], [AS_IF([test ! -x "$with_bowtie2" || test -d "$with_bowtie2"], [AC_MSG_ERROR([Error: bowtie2 path does not point to an executable.])], [bowtie2check=$with_bowtie2])])
AC_ARG_ENABLE([bowtie2], AS_HELP_STRING([--disable-bowtie2], [Disable bowtie2 usage, NanoVar will not work with short-reads.]))
AS_IF([test "x$enable_bowtie2" = xno], [echo Disabled bowtie2, NanoVar will not work with short-reads.], [AS_IF([test "x$bowtie2check" = xno], [AC_MSG_ERROR([Error: bowtie2 not found, please ensure "bowtie2" in \$PATH before installing NanoVar, or ./configure with --with-bowtie2=/path/to/bowtie2, or ./configure with --disable-bowtie2])])])
AS_IF([test "x$bowtie2check" != xno], [bowversion=$($bowtie2check --version | sed 's|.*version ||g' | head -n 1)])
AC_SUBST([bowversion])

#bowtie2-build check
AC_ARG_WITH([bowtie2build], AS_HELP_STRING([--with-bowtie2build=/path/to/bowtie2-build], [Specify full path to bowtie2-build]), [echo Assigning bowtie2-build... "$with_bowtie2build"], [AC_PATH_PROG(bowtie2buildcheck,bowtie2-build,no)])
AS_IF([test ! -z "$with_bowtie2build"], [AS_IF([test ! -x "$with_bowtie2build" || test -d "$with_bowtie2build"], [AC_MSG_ERROR([Error: bowtie2-build path does not point to an executable.])], [bowtie2buildcheck=$with_bowtie2build])])
AS_IF([test "x$enable_bowtie2" = xno], [echo Disabled bowtie2-build, NanoVar will not work with short-reads.], [AS_IF([test "x$bowtie2buildcheck" = xno], [AC_MSG_ERROR([Error: bowtie2-build not found, please ensure "bowtie2-build" in \$PATH before installing NanoVar, or ./configure with --with-bowtie2build=/path/to/bowtie2-build, or ./configure with --disable-bowtie2])])])

AM_CONDITIONAL([DISABLE_BOWTIE2], [test "x$enable_bowtie2" = xno])

AC_CONFIG_FILES([Makefile scripts/Makefile])
AC_OUTPUT
