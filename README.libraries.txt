README.libraries.txt:
=====================

You maybe need to compile the code against some external
libraries. This lately happened in the case of parton distribution
functions in the context of PYTHIA calculations. 
The following explains the chosen concept using this example. It is
straightforward to extend it to other cases.

For the impatient:
------------------
You may compile the code with "make PDF=PDFLIB" or "make PDF=LHAPDF"
or just "make". The latter just uses some stub routines and is the old
behavior. If you want to use one of the two first methods, you have
to read the whole scrap.

For the more patient:
---------------------
As known all *.o files are kept in the directory /objects. We include
now the directory /objects/LIB/lib in the compile and link command via
"-L/objects/LIB/lib". (The doubling of "lib" was necessary, since some
'configure; make; make install' chains for some libraries create
additional directories, which end for example as /objects/LIB/bin,
/objects/LIB/include etc.) For the special case considered here, we
also include "-lPDF" as argument for the compiler. This means, we
have to provide /objects/LIB/lib/libPDF.a, which we choose to be a
soft link to the actual library to use. These are:
* /objects/LIB/lib/libLHAPDF.a
* /objects/LIB/lib/libPDF8.04.a 
* /objects/LIB/lib/libPDFstub.orig.6225.a

The first two libraries you have to provide (if desired) by your own,
the third one is generated automatically by the make procedure.

How is libPDFstub.orig.6225.a generated?
----------------------------------------
This library is converted from libPDFstub.orig.6225.o, which is an
ordinary object file from the compilation process. We prohibit the *.o
file to be included in object file lists and to be used in the
linking. This library remains in /objects and only a soft link in
/objects/LIB/lib is created.


(--- THE FOLLOWING IS ONLY FOR NON-PUBLIC-RELEASE-USERS ---)


How to provide external libraries
---------------------------------
We have prepared a directory tree /libraries, which runs parallel to
/workingCode. Here you can add the tarballs etc from the used
libraries and compile them as usual. (We will not provide any of them
in the public releases of the GiBUU code.) For the actual case we have
* libraries/PDF/lhapdf-5.4.0
* libraries/PDF/pdflib-8.04
as directories keeping the provided libraries.

LHAPDF is a library with a full 'configure; make; make install'
chain. PDFLIB8.04 is an own cut-out off the CERNLIB with a simple
Makefile structure. 
In /libraries/PDF you find a Makefile, which implements the necessary
calls for the libraries. BUT: in order to build the libraries with the
necessary compiler selection and flags, you have to initiate it in the
directory /workingCode. (We have provided some Makefile targets to do
this for you. We also avoid the pitfall, that the "working code"
directory is a soft link to some \workingCode somewhere.)