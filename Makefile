# Makefile
# This file is part of LIME, the versatile line modeling engine
#
# Copyright (C) 2006-2014 Christian Brinch
# Copyright (C) 2015-2017 The LIME development team

# Platform-dependent stuff:
include Makefile.defs

##
## Make sure to put the correct paths.
##
PREFIX  =  ${PATHTOLIME}

# Paths:
srcdir		= ${CURDIR}/src
docdir		= ${CURDIR}/doc
exampledir	= ${CURDIR}/example
#*** better to use ${PREFIX} here rather than ${CURDIR}? (the latter is used in artist/lime.)

ifneq (,$(wildcard ${PREFIX}/lib/.))
    LIBS += -L${PREFIX}/lib
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


CPPFLAGS += -I${PREFIX}/include \
	    -I${PREFIX}/src \
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
PYTARGET = pylime
CC	= gcc -fopenmp
MODELS  = model.c # Overwritten in usual practice by the value passed in by the 'lime' script.
MODELO 	= ${srcdir}/model.o

CCFLAGS = -O3 -falign-loops=16 -fno-strict-aliasing
LDFLAGS = -lgsl -lgslcblas -l${LIB_QHULL} -lcfitsio -lncurses -lm 

ifeq (${DOTEST},yes)
  CCFLAGS += -DTEST
#  CC += -g -Wunused -Wno-unused-value -Wformat -Wformat-security
  CC += -g -Wall
endif

ifeq (${VERBOSE},no)
  CCFLAGS += -DNOVERBOSE
endif

ifeq (${USEHDF5},yes)
  CPPFLAGS += -DUSEHDF5
  CCFLAGS += -DH5_NO_DEPRECATED_SYMBOLS
  LDFLAGS += -lhdf5_hl -lhdf5 -lz
  CORESOURCES += ${HDF5SOURCES}
  CONVSOURCES += ${HDF5SOURCES}
  COREINCLUDES += ${HDF5INCLUDES}
endif

SRCS = ${CORESOURCES} ${STDSOURCES}
INCS = ${COREINCLUDES}
PYSRCS = ${CORESOURCES} ${PYSOURCES}
PYINCS = ${COREINCLUDES} ${PYINCLUDES}
OBJS = $(SRCS:.c=.o)
PYOBJS = $(PYSRCS:.c=.o)
CONV_OBJS = $(CONVSOURCES:.c=.o)

.PHONY: all doc docclean objclean limeclean clean distclean pyclean python

all:: ${TARGET} 

# Implicit rules:
%.o : %.c
	${CC} ${CCFLAGS} ${CPPFLAGS} -o $@ -c $<

${TARGET}: ${OBJS} ${MODELO} 
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}

${OBJS} : ${INCS}
${CONV_OBJS} : ${CONVINCLUDES}
${PYOBJS} : ${PYINCS}

${MODELO}: ${INCS}
	${CC} ${CCFLAGS} ${CPPFLAGS} -o ${MODELO} -c ${MODELS}

python: CCFLAGS += ${PYCCFLAGS}
python: CPPFLAGS += -DNO_NCURSES
python: LDFLAGS += ${PYLDFLAGS}
python: ${PYTARGET}

${PYTARGET}: ${PYOBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}

gridconvert : CPPFLAGS += -DNO_NCURSES

gridconvert: ${CONV_OBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}

doc::
	mkdir ${docdir}/_html || true
	sphinx-build doc ${docdir}/_html

docclean::
	rm -rf ${docdir}/_html

objclean::
	rm -f *~ ${srcdir}/*.o

limeclean:: objclean
	rm -f ${TARGET}

pyclean:: objclean
	rm -f ${pydir}/*.pyc ${PYTARGET}

clean:: objclean pyclean
	rm -f gridconvert

distclean:: clean docclean limeclean
	rm Makefile.defs

