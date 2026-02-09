#
# A simple makefile that links Absoft Fortran 90 against the IMSL library
#

FC = f90
FOPTS = -O3 -cpu:p7 -N113 -YEXT_NAMES=LCS -YEXT_SFX=_

OUT = code24.e
SRC = code24.f

$(OUT) : code24.f code24.o
	$(FC) -o $(OUT) code24.o -L/usr/absoft/lib -limsl -limslblas

.f.o:
	$(FC) -c $(FOPTS) $*.f
