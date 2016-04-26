# Makefile for convert
FC = gfortran

all: remstart compile clean

compile:
	$(FC) convert.f -o convert.x
	$(FC) grid_interpol.f -o grid_interpol.x

remstart:
	rm -rf *.o convert.x grid_interpol.x
	
clean:
	rm -rf *.o
