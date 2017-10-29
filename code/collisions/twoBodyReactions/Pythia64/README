                  ****************************
                  *                          *
                  *          README          *
                  *                          *
                  *    PYTHIA version 6.4    *
                  *                          *
                  ****************************

                  (Last updated 03 April 2007)

PYTHIA version 6.4 is a direct continuation of version 6.3;
actually 6.400 is identical with 6.327. Therefore it should not
be a big operation for the normal user to run the program.

INSTALLATION NOTES:
   1. Previously, PYTHIA was distributed as one large file.  This
format is not available from svn, but is from the downloads area
of the Pythia6 HepForge page.  

   2. In the current distribution, it is expected that the
user will make a library that is then linked with their own
main program and analysis code.  

   3. We recommend that the PYTHIA tarball be unpacked in
a separate directory, e.g. "mkdir pylib; tar xzf pythia-6.4.10.tgz -C pylib"

   4. To make the library, type "make" for options.
       a. "make clean" will remove the library and the object files
       b. "make lib" will make the library
       c. Rather than typing out all of the subroutine and function
          names, which may grow with future distributions and otherwise
          break things, the Makefile uses wildcard substitution, so it is 
          important that you have no private code in the same directory 
          that conflicts with the names in the Makefile, e.g. "py*.f ..."

   5. After "make lib", a file "libpythia.a" should exist.

   6. When building your main code, you should then link in this library,
e.g. "f77 -c main.exe main.f -L$(path to pylib) -lpythia"

   7. In the main program, the fortran command "EXTERNAL PYDATA" will ensure
that the common block default values can be changed by the user. (*)

   8. This Makefile has been tested on a Linux distribution only.


(*) An EXTERNAL statement is used to identify a symbolic name as representing an 
external procedure or dummy procedure, and to permit such a name to be used as an 
actual argument.  The form of an EXTERNAL statement is:

EXTERNAL proc [,proc]...

where each proc is the name of an external procedure, dummy procedure, or block data 
subprogram. Appearance of a name in an EXTERNAL statement declares that name to be an 
external procedure name or dummy procedure name, or block data subprogram name. 

   -- from www.fortran.com/fortran/F77_std
