BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

if FRAME
if LALFRAME
FRAMESRCS = \
	series.c \
	series.h \
	lalappsfrutils.c \
	lalappsfrutils.h
endif
endif

if LALMETAIO
LALMETAIOSRCS = \
	processtable.c \
	processtable.h
endif

if FFTW
lalapps_fftwf_wisdom_SOURCES = \
	fftwf_wisdom.c
lalapps_fftw_wisdom_SOURCES = \
	fftw_wisdom.c
FFTWPROGS = \
	lalapps_fftwf_wisdom \
	lalapps_fftw_wisdom
endif

noinst_LTLIBRARIES = liblalapps.la
liblalapps_la_SOURCES = \
	lalapps.c \
	lalapps.h \
	LALAppsVCSInfo.c \
	LALAppsVCSInfo.h \
	$(FRAMESRCS) \
	$(LALMETAIOSRCS)
nodist_liblalapps_la_SOURCES = LALAppsBuildInfo.c

noinst_HEADERS = \
	LALAppsVCSInfoHeader.h

bin_PROGRAMS = \
	lalapps_version \
	lalapps_tconvert \
	lalapps_cache \
	$(FFTWPROGS)

lalapps_version_SOURCES = \
	version.c \
	lalapps.h \
	LALAppsVCSInfo.h
lalapps_tconvert_SOURCES = \
	tconvert.c \
	getdate.c \
	getdate.h
lalapps_cache_SOURCES = \
	cache.c

LDADD = liblalapps.la

if HAVE_PYTHON
pkgpython_PYTHON = \
	__init__.py \
	git_version.py
pybin_scripts = \
	lalapps_path2cache \
	lalapps_searchsum2cache
endif

TESTS = lalapps_version

vcs_info_sources = LALAppsVCSInfoHeader.h git_version.py
build_info_source = LALAppsBuildInfo.c
$(build_info_source): $(liblalapps_la_DEPENDENCIES)

EXTRA_DIST += \
	__init__.py \
	fftw_wisdom.c \
	fftwf_wisdom.c \
	getdate.y \
	git_version.py \
	lalapps.dox \
	lalappsfrutils.c \
	lalappsfrutils.h \
	processtable.c \
	processtable.h \
	series.c \
	series.h \
	$(END_OF_LIST)
