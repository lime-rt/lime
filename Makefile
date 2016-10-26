# Makefile
# This file is part of LIME, the versatile line modeling engine
#
# Copyright (C) 2006-2014 Christian Brinch
# Copyright (C) 2015-2016 The LIME development team

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

CPPFLAGS	= -I${PREFIX}/include \
		  -I${PREFIX}/src \
		  -I${HOME}/include \
		  -I/opt/local/include \
		  -I/sw//include \
	          ${EXTRACPPFLAGS}

ifdef OLD_QHULL
	QHULL   = qhull
	CPPFLAGS += -DOLD_QHULL
else
	QHULL   = qhullstatic
endif

ifdef OLD_FITSIO
	CPPFLAGS += -DOLD_FITSIO
endif

# Names of source files included:
include Makefile.defs

##
## Do not change anything below unless you know what you are doing! 
##

TARGET  = lime.x # Overwritten in usual practice by the value passed in by the 'lime' script.
PYTARGET = pylime
CC	= gcc -fopenmp
MODELS  = model.c # Overwritten in usual practice by the value passed in by the 'lime' script.
MODELO 	= src/model.o

CCFLAGS = -O3 -falign-loops=16 -fno-strict-aliasing
LDFLAGS = -lgsl -lgslcblas -l${QHULL} -lcfitsio -lncurses -lm 

#vvvvv
# These will be installation-dependent and should therefore be set up via autoconf. The flag values can be obtained for any given system by running
#
#	<python prefix>-config --cflags
# and
#	<python prefix>-config --ldflags
#
# respectively.
PYCCFLAGS = -I/usr/include/python2.7 -I/usr/include/python2.7 -fno-strict-aliasing -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 -m64 -mtune=generic -D_GNU_SOURCE -fPIC -fwrapv -DNDEBUG -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 -m64 -mtune=generic -D_GNU_SOURCE -fPIC -fwrapv

PYLDFLAGS = -lpython2.7 -lpthread -ldl -lutil -lm -Xlinker -export-dynamic
#^^^^^

ifeq (${DOTEST},yes)
  CCFLAGS += -DTEST
  CC += -g
endif

SRCS = ${CORESOURCES} ${STDSOURCES}
INCS = ${COREINCLUDES}
OBJS = $(SRCS:.c=.o)
PYSRCS = ${CORESOURCES} ${PYSOURCES}
PYINCS = ${COREINCLUDES} ${PYINCLUDES}
PYOBJS = $(PYSRCS:.c=.o)

.PHONY: all clean distclean python
	all:: ${TARGET} 

# Implicit rules:
%.o : %.c
	${CC} ${CCFLAGS} ${CPPFLAGS} -o $@ -c $<

${TARGET}: ${OBJS} ${MODELO} 
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}

${OBJS} : ${SRCS} ${INCS}

${MODELO}: ${MODELS} ${INCS}
	${CC} ${CCFLAGS} ${CPPFLAGS} -o ${MODELO} -c ${MODELS}

python: CCFLAGS += ${PYCCFLAGS}
python: CPPFLAGS += -DNO_NCURSES
python: LDFLAGS += ${PYLDFLAGS}
python: ${PYTARGET}

${PYTARGET}: ${PYOBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}

${PYOBJS} : ${PYSRCS} ${PYINCS}

doc::
	mkdir ${docdir}/_html || true
	sphinx-build doc ${docdir}/_html

docclean::
	rm -rf ${docdir}/_html

clean:: 
	rm -f *~ ${srcdir}/*.o ${pydir}/*.pyc ${TARGET} ${PYTARGET}

distclean:: clean docclean

