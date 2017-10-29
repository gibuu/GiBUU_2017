.SUFFIXES:            # Delete the default suffixes

all: ReadArr.so

tmp.f90: Read.F90
	cpp -o tmp.f90 -E Read.F90

ReadArr.so: tmp.f90
	f2py -c -m ReadArr tmp.f90