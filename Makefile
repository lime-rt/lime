##
## Make sure to put the correct paths.
##
PREFIX  =  ${PATHTOLIME}
LIBS	= -L${PREFIX}/lib \
		  -L${HOME}/lib \
          -L/opt/local/lib 
INCLUDE	= -I${PREFIX}/include \
		  -I${PREFIX}/src \
		  -I${HOME}/include \
		  -I/opt/local/include

# For qhull version 2011.1 or newer:
QHULL   = qhullstatic
# For older qhull versions:
#QHULL   = qhull



##
## Do not change anything below unless you know what you are doing! 
##

TARGET  = lime.x 
CC		= gcc
SRCS    = src/aux.c src/curses.c src/grid.c src/LTEsolution.c   \
		  src/main.c src/molinit.c src/photon.c src/popsout.c   \
		  src/predefgrid.c src/ratranInput.c src/raytrace.c     \
		  src/smooth.c src/sourcefunc.c src/stateq.c            \
          src/statistics.c src/magfieldfit.c src/stokesangles.c \
		  src/writefits.c src/weights.c src/velospline.c        \
		  src/old_raytrace.c src/getclosest.c
MODELS  = /Users/christianbrinch/Work_projects/Jes/irs43+63/model/lime/irs63.c
OBJS    = src/aux.o src/curses.o src/grid.o src/LTEsolution.o   \
		  src/main.o src/molinit.o src/photon.o src/popsout.o   \
		  src/predefgrid.o src/raytrace.o src/ratranInput.o     \
		  src/smooth.o src/sourcefunc.o src/stateq.o            \
		  src/statistics.o src/magfieldfit.o src/stokesangles.o \
		  src/writefits.o src/weights.o src/velospline.o        \
		  src/old_raytrace.o src/getclosest.o
MODELO 	= src/model.o

CCFLAGS = -O3 -falign-loops=16 -fno-strict-aliasing   
LDFLAGS = -lgsl -lgslcblas -l${QHULL} -lcfitsio -lncurses -lm 

.SILENT:

.PHONY: all clean distclean 
	all:: ${TARGET} 

${TARGET}: ${OBJS} ${MODELO} 
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}   

${MODELO}:  
	${CC} ${CCFLAGS} ${INCLUDE} -o ${MODELO} -c ${MODELS} 

${OBJS}: %.o: %.c  
	${CC} ${CCFLAGS} ${INCLUDE} -o $@ -c $< 

clean:: 
	rm -f *~ src/*.o ${TARGET} 

distclean:: clean

