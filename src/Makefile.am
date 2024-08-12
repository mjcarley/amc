AM_CPPFLAGS = -I$(top_srcdir) -I$(includedir) -DG_LOG_DOMAIN=\"AMC\"

lib_LTLIBRARIES = libamc.la

libamc_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

libamc_la_SOURCES = \
	amc.c \
	amc.h \
	frame.c \
	tinyexpr.c \
	tinyexpr.h

include_HEADERS = \
	amc.h

CLEANFILES = $(BUILT_SOURCES)
