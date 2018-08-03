# Makefile
# This file is part of LIME, the versatile line modeling engine
#
# Copyright (C) 2006-2014 Christian Brinch
# Copyright (C) 2015-2017 The LIME development team

# Platform-dependent stuff:
include Makefile.defs

# Paths:
srcdir		= ${CURDIR}/src
docdir		= ${CURDIR}/doc
exampledir	= ${CURDIR}/example
pydir		= ${CURDIR}/python

ifneq (,$(wildcard ${CURDIR}/lib/.))
    LIBS += -L${CURDIR}/lib
endif
ifneq (,$(wildcard ${HOME}/lib/.))
    LIBS += -L${HOME}/lib
endif
ifneq (,$(wildcard /opt/local/lib/.))
    LIBS += -L/opt/local/lib
endif
ifneq (,$(wildcard /sw/lib/.))
    LIBS += -L/sw/lib
endif
ifneq (,$(wildcard /usr/local/lib/.))
    LIBS += -L/usr/local/lib
endif


CPPFLAGS += -I${CURDIR}/include \
	    -I${CURDIR}/src \
	    -I${HOME}/include \
	    -I/opt/local/include \
	    -I/sw//include \
	    ${EXTRACPPFLAGS}


# Names of source files included:
include Makefile.srcs

##
## Do not change anything below unless you know what you are doing! 
##

TARGET  = lime.x # Overwritten in usual practice by the value passed in by the 'lime' script.
MLRUMP = modellib
LLRUMP = lime
MLTARGET = lib${MLRUMP}.so
LLTARGET = lib${LLRUMP}.so
CC	= gcc -fopenmp
MODELO 	= ${srcdir}/model.o

CCFLAGS += -O3 -falign-loops=16 -fno-strict-aliasing
LDFLAGS += -lgsl -lgslcblas -l${LIB_QHULL} -lcfitsio -lm 

ifeq (${DOTEST},yes)
  CCFLAGS += -DTEST
#  CC += -g -Wunused -Wno-unused-value -Wformat -Wformat-security
  CC += -g -Wall
endif

ifeq (${VERBOSE},no)
  CCFLAGS += -DNOVERBOSE
endif

ifeq (${USECURSES},yes)
  LDFLAGS += -lncurses
  CCFLAGS += -DNCURSES
endif

ifeq (${USEHDF5},yes)
  CPPFLAGS += -DUSEHDF5
  CCFLAGS += -DH5_NO_DEPRECATED_SYMBOLS
  LDFLAGS += -lhdf5_hl -lhdf5 -lz
  CORESOURCES += ${HDF5SOURCES}
  CONVSOURCES += ${HDF5SOURCES}
  COREINCLUDES += ${HDF5INCLUDES}
else
  CORESOURCES += ${FITSSOURCES}
  CONVSOURCES += ${FITSSOURCES}
  COREINCLUDES += ${FITSINCLUDES}
endif

SRCS = ${CORESOURCES} ${STDSOURCES}
INCS = ${COREINCLUDES}
OBJS = $(SRCS:.c=.o)

ifeq (${USECURSES},yes)
  SRCS += ${CUR_SOURCE}
else
  SRCS += ${MSG_SOURCE}
endif

PYSRCS = ${CORESOURCES} ${PYSOURCES}
PYINCS = ${COREINCLUDES} ${PYINCLUDES}
PYOBJS = $(PYSRCS:.c=.o)

MLSRCS = ${MLCORESOURCES} ${MLSOURCES}
MLINCS = ${COREINCLUDES} ${MLINCLUDES}
MLOBJS = $(MLSRCS:.c=.o)

LLSRCS = ${CORESOURCES} ${LLSOURCES}
LLINCS = ${COREINCLUDES}
LLOBJS = $(LLSRCS:.c=.o)

CLSRCS = ${CORESOURCES} ${MLCORESOURCES} ${CLSOURCES}
CLINCS = ${COREINCLUDES}
CLOBJS = $(CLSRCS:.c=.o)

CONV_OBJS = $(CONVSOURCES:.c=.o)

.PHONY: all doc docclean objclean limeclean clean distclean pyclean pyshared pyshclean casa casaclean

all:: ${TARGET} 

# Implicit rules:
%.o : %.c
	${CC} ${CCFLAGS} ${CPPFLAGS} -o $@ -c $<

${OBJS} : ${INCS}
${PYOBJS} : ${PYINCS}
${MLOBJS} : ${MLINCS}
${LLOBJS} : ${LLINCS}
${CLOBJS} : ${CLINCS}
${CONV_OBJS} : ${COREINCLUDES}

${MODELO}: ${INCS}
	${CC} ${CCFLAGS} ${CPPFLAGS} -o ${MODELO} -c ${MODELS}

${TARGET}: ${MODELO} ${OBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}

pylime: CCFLAGS += ${PYCCFLAGS}
pylime: CPPFLAGS += -DIS_PYTHON
pylime: LDFLAGS += ${PYLDFLAGS}

pylime: ${PYOBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}
	rm -f ${srcdir}/*.o ${srcdir}/*/*.o

casalime: CCFLAGS  += ${PYCCFLAGS}
casalime: CPPFLAGS += -DIS_PYTHON -DNO_PROGBARS
casalime: LDFLAGS  += ${PYLDFLAGS}

casalime: ${CLOBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}
	rm -f ${srcdir}/*.o ${srcdir}/*/*.o

${MLTARGET}: ${MLOBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}

${LLTARGET}: LIBS += -L${CURDIR}

# This way liblime.so can always find libmodellib.so without the user needing to set LD_LIBRARY_PATH. The extra $ seems to be an escape character needed by make; the command output string should be
#	-Wl,-rpath,'$ORIGIN'
#
XLDFLAGS = -Wl,-rpath,'$$ORIGIN'
${LLTARGET}: LDFLAGS += ${XLDFLAGS}

${LLTARGET}: ${MLTARGET} ${LLOBJS}
	${CC} -o $@ ${LLOBJS} ${LIBS} ${LDFLAGS} -l${MLRUMP}

pyshared: CCFLAGS  += ${PYCCFLAGS} -fPIC
pyshared: CPPFLAGS += -DIS_PYTHON
pyshared: LDFLAGS  += ${PYLDFLAGS} -shared

pyshared: ${LLTARGET}

gridconvert: ${CONV_OBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}

doc::
	mkdir ${docdir}/_html || true
	sphinx-build doc ${docdir}/_html

docclean::
	rm -rf ${docdir}/_html

objclean::
	rm -f ${srcdir}/*.o ${srcdir}/*/*.o

limeclean:: objclean
	rm -f ${TARGET}

pyclean:: objclean
	rm -f ${pydir}/*.pyc pylime

pyshclean:: objclean
	rm -f ${MLTARGET} ${LLTARGET}

casaclean:: objclean
	rm -f casalime

clean:: objclean pyclean pyshclean casaclean
	rm -f gridconvert
	rm -f *~ ${srcdir}/*~ ${srcdir}/*/*~

distclean:: clean docclean limeclean
	rm Makefile.defs

