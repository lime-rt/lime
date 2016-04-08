# Makefile
# This file is part of LIME, the versatile line modeling engine
#
# Copyright (C) 2006-2014 Christian Brinch
# Copyright (C) 2015 The LIME development team

##
## Make sure to put the correct paths.
##
PREFIX  =  ${PATHTOLIME}

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

##
## Do not change anything below unless you know what you are doing! 
##

TARGET  = lime.x 
CC		= gcc -fopenmp
SRCS    = src/aux.c src/messages.c src/grid.c src/LTEsolution.c   \
		  src/main.c src/molinit.c src/photon.c src/popsin.c    \
		  src/popsout.c src/predefgrid.c src/ratranInput.c      \
          src/raytrace.c src/smooth.c src/sourcefunc.c          \
		  src/stateq.c src/statistics.c src/magfieldfit.c       \
		  src/stokesangles.c src/writefits.c src/weights.c      \
		  src/velospline.c src/getclosest.c src/raythrucells.c  \
		  src/tcpsocket.c src/defaults.c src/fastexp.c
MODELS  = model.c
OBJS    = src/aux.o src/messages.o src/grid.o src/LTEsolution.o   \
		  src/main.o src/molinit.o src/photon.o src/popsin.o    \
		  src/popsout.o src/predefgrid.o src/raytrace.o         \
		  src/ratranInput.o src/smooth.o src/sourcefunc.o       \
		  src/stateq.o src/statistics.o src/magfieldfit.o       \
		  src/stokesangles.o src/writefits.o src/weights.o      \
		  src/velospline.o src/getclosest.o src/raythrucells.o  \
		  src/tcpsocket.o src/defaults.o src/fastexp.o
MODELO 	= src/model.o

#CCFLAGS = -O3 -falign-loops=16 -fno-strict-aliasing -DTEST -DFASTEXP -DNO_NCURSES
CCFLAGS = -O3 -falign-loops=16 -fno-strict-aliasing
LDFLAGS = -lgsl -lgslcblas -l${QHULL} -lcfitsio -lncurses -lm 

.SILENT:

.PHONY: all clean distclean 
	all:: ${TARGET} 

${TARGET}: ${OBJS} ${MODELO} 
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}  

${MODELO}:  
	${CC} ${CCFLAGS} ${CPPFLAGS} -o ${MODELO} -c ${MODELS}

${OBJS}: %.o: %.c  
	${CC} ${CCFLAGS} ${CPPFLAGS} -o $@ -c $<

clean:: 
	rm -f *~ src/*.o ${TARGET} 

distclean:: clean

