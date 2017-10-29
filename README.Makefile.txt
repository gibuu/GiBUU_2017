#****e* /README.Makefile.txt
# NAME
# README.Makefile.txt
#
# SOURCE
#

This is a README file concerning Makefile... 

[...This text is under construction. Please modify! ...]






FOR THE IMPATIENT:
==================

Just type
>>  make
in order to create all "*.o" files and and an executable "GiBUU.x",
which is located in "testrun".

(If this does not work, while the complaints are about missing "make
targets", type
>>  make renew
and retry.)

To run the program, just type
>>  cd testrun
>>  ./GiBUU.x < "Name of your jobcard"

If you want to delete all compilation output, type
>>  make veryclean

If you want to delete all compilation output and all kind of other 
spurious files (*~,*.dat,#*,fort.*) , type
>>  make superclean


DETAILLED INFORMATION:
======================

1) The Philosophy of (GNU) make:
--------------------------------
The GNU make command allows with the aid of a "Makefile" to build all
necessary files in the right order, i.e. to respect all dependencies
between the different source code files.
Unfortunately, (GNU) make is not able to hande dependencies between
the source code files of different directories. Especially, it is not
capable to handle such a detailed directory tree structure as invented
her for the BUU code.

2) The Philosophy of the BUU Makefile:
--------------------------------------
Because (GNU) make can not handle dependencies between different
directories automatically, the programmers have to take care for it
themselves. 

We have choosen the philosophy, that we provide in the main directory
a file "Makefile" and some templates. Typing "make" calls the
"Makefile", which then first (by calling the subtarget "make renew")
copies some Makefiles everywhere in the tree structure, where it is
necessary. Then in a second step (within the same main call), the make
command iterates over all subdirectories in the tree and creates
(soft-)links of the source files to some unique directory (i.e
"objects/"). The third step is to create the dependencies between the
different source files. This is either done by "makedepf90" or an own
Perl-Script called "Own_Makedepf90.pl". And finally, in the fourth
step, the compiler is called to compile all source files and link them
to an executable called "GiBUU.x"

3) The Makefile templates and where to copy them:
-------------------------------------------------
While iterating over the source tree, the major "make" has to know,
which subdirectories it should step in.
In order to catch the complex structures, you have to provide in every
subdirectory a file named "Makefile.local". Whith this you get some
influence of the traversal of the tree:

1) Normally, make searches subdirectories in alphabetical order, 
BUT:
1a) subdirectories given via the variable "SUBDIR" are searched first
    and in the given order.
1b) subdirectories given via the variable "SUBDIREXCL" are excluded
    from the directory list.
1c) subdirectories not met by these to variables are included
    automatically in alphabetically ordering.

2) Normally, make compiles all ".f" and "*.f90" files, 
BUT:
2a) filenames given via the variable "FILEEXCL" are excluded from the
    compile list.



4) How to use:
--------------

0) If your source code tree does not provide any Makefile in the
subdirectories, or if you want to replace all Makefiles in the tree, 
please type
>> make renew

This command will delete all your *.o and *.mod files, but also
generate a corresponding Makefile in every subdirectory.

Please note: 
typing
>> make
tries to find out, whether all the Makefiles have to be rebuild or not.

Any error messages at this stage show already problems you have to cure.

1) in order to produce all files, type
>>  make

Error messages at this stage are mainly due to 
- programming errors
- faulty dependencies!



5) What you can get from the Makefile:
--------------------------------------

The main Makefile knows the following targets:

-) make
   recompile all necessary *.o and *.mod and produce "main.x"

-) make clean
   delete: *.~

-> make cleanEXE
   a) delete *.x in every subdirectory
   b) delete *.o and *.mod in code/objects
   (forces a new compilation of main.x. *.o/*.mod remain in the
    corresponing subdirectories, no recompilation necessary.) 

-) make veryclean
   delete: *.~ *.o *.mod *.x in every subdirectory
   (forces a new compilation of all files)

-) make MAKEFILEclean:
   delete: *.~ *.o *.mod *.x. MakefileDepend

-) make Makefiles
   copy all templates to all the subdirectories

-) make renew:
   abbrev. for:  make MAKEFILEclean; make MAKEFILE

-) make version
   write svn version information of the directory to file 
   version.f90 such that code is returning those infos

-) make quick
   = "make" without "make makeVersion"

=> More: There are event more Makefile entries available: 
   For an up-to-date list of all entries see the HTML-Documentation!


7) Alternative Compilers:
-------------------------

By default the Intel Compiler ("ifort") is used to compile GiBUU,
but also the GCC/gfortran, g95, Sun/Oracle and Portland compilers
are supported. To use one of these you have to set the FORT variable:

make FORT=gfortran

This variable should contain the name of the compiler executable
which you want to use (and which should be present in your PATH).
Possible choices are e.g.

FORT={ifort,gfortran,g95,sunf95,pgf95, ...}


8) Additional Options:
----------------------

By default the GiBUU code is compiled with debugging flags.
If you don't need debugging info or want an optimized executable,
you can set the MODE variable:

make MODE=opt

This compiles the code with the "-O3" flag (highest optimization).
In the same way you can use MODE={opt1,opt2,opt3} to control the
optimization level (-O1, -O2, -O3).

In addition to optimization, you can also tell the compiler to
produce profiling information, which can be used with a profiler
like "gprof". This is done via:

make MODE=prof

For the Intel compiler also the level of floating point exceptions
can be adjusted. The default is -fpe3. This can be changed with e.g.

make FPE=0

Possible values here are FPE={0,1,2,3}. See Intel compiler
documentation for details.

9) parton distribution libraries
--------------------------------

You may switch to link against some parton distribution libraries (if
provided). This is done via

make PDF=PDFLIB

or

make PDF=LHAPDF

if you want to switch to the corresponding libaray (see
README.libraries.txt for more information).

#*******************
