#
# A simple makefile that links Absoft Fortran 90 against the IMSL library
#

FC = f90
FOPTS = -O3 -cpu:p7 -N113 -YEXT_NAMES=LCS -YEXT_SFX=_

OUT = tmp.e
SRC = tmp.f

$(OUT) : tmp.f tmp.o
	$(FC) -o $(OUT) tmp.o -L/usr/absoft/lib -limsl -limslblas

.f.o:
	$(FC) -c $(FOPTS) $*.f
