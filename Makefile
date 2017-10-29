#*************************************************** -*- Makefile -*- **********
#****e* /Makefile
# NAME
# Makefile
# PURPOSE
# This is the main Makefile of the GiBUU code.
#
# you may compile your code via
#   make [VAR=val | ...]
# with VAR = ...
#
# FORT:
# * FORT = ifort,gfortran,...
# * FORT = /path/to/intel/compiler/ifort
#
# MODE:
# * MODE = opt0,opt1,opt2,opt3
# * MODE = opt4
# * MODE = opt5
# * MODE = lto
# * MODE = prof
# * MODE = callGraph
#
# STATIC:
# * STATIC = 0|1
#
# FPE:
# * FPE = 0|1|2|3
#
# ARCH:
# * ARCH = 32
#
# ARGS:
# * ARGS = "..."
#
# NOTES
# The most aggresive optimization you can get with:
# * make FORT=gfortran MODE=lto ARGS="-fopenmp -march=native"
# * make FORT=ifort MODE=opt3 ARGS="-fopenmp -parallel -xHost"
#*******************************************************************************





#*******************************************************************************
#*******************************************************************************
#************* 1) Setting variables
#*******************************************************************************
#*******************************************************************************

.SUFFIXES:            # Delete the default suffixes

export SHELL = /bin/bash

### Name of GiBUU main file
export NameOfExe=GiBUU

export OS = $(strip $(shell uname))
export OS_LONG = '$(shell uname -s -r -m)'

export noPrintDirectory=--no-print-directory

### PROGRAMS:
export ECHO = echo -e
export PERL = $(strip $(shell which perl 2>/dev/null))
ifeq ($(PERL),$(EMPTY))
  export MAKEDEP = $(strip $(shell which makedepf90 2>/dev/null))
endif
ifeq ($(OS),Darwin)
  export FIND = gfind
else
  export FIND = find
endif
export AR = ar

### DIRECTORIES:

export ROOTDIR := $(CURDIR)

# if CURDIR is a symbolc link, then use better the following:
export ROOTDIRLONG := `pwd -P`

# Directory to store the *.o and *.mod files:
export OBJDIR := $(CURDIR)/objects

# Directory to store the executable:
export EXEDIR := $(CURDIR)/testRun

# Directory to find sources of additional libraries:
export LIBDIRSRC := $(ROOTDIRLONG)/../libraries

# Directory to find the .a files:
export LIBDIR := $(OBJDIR)/LIB/lib

### Input directory

export PATH_TO_INPUT='../buuinput'

# The variable "allSrcFiles" includes :
# * path of all the source files which are linked to objects
# * path to the main file.
SRCf90 := $(wildcard objects/*.f90)
SRCF90 := $(wildcard objects/*.F90)
SRCf77 := $(wildcard objects/*.f)
allSrcFiles:= $(SRCf90) $(SRCF90) $(SRCf77) code/$(NameOfExe).f90
export allSrcFiles:= $(filter-out objects/version.f90, $(allSrcFiles))

# Format definition for "nicer" outputs #######
# Use e.g. :
#   "make StartHeader= EndHeader= StartHeader_green= StartHeader_red= blue= endBlue= "
# to override the formats. Possible colors for bash are:
#export black='\033[30;47m'
#export red='\033[31;47m'
redWhite='\033[31m'
#export green='\033[32;47m'
whiteGreen='\033[37;42m'
whiteRed='\033[37;41m'
#export greenWhite='\033[32m'
#export yellow='\033[33;47m'
#export yellowWhite='\033[33m'
#export blueWhite='\033[34m'
#export magenta='\033[35;47m'
#export cyan='\033[36;47m'
#export white='\033[37;47m'
export blue='\033[34;47m'
export endBlue='\033[0m'
export StartHeader=$(redWhite)"\033[1m
export StartHeader_green=$(whiteGreen)"\033[1m
export StartHeader_red=$(whiteRed)"\033[1m
export EndHeader=\033[0m"



### SUBDIRS:

SUBDIR := code
#*******************************************************************************
#*******************************************************************************
#************* 2) Defining the compiler
#*******************************************************************************
#*******************************************************************************

EMPTY =

### do dynamic linking by default
STATIC = 0

### build with optimizations by default
MODE = opt

FPE = 3

### default compiler: Intel
FORT=ifort
FORTPATH = $(strip $(shell which $(FORT) 2>/dev/null))

### fallback options if ifort is not available
ifeq ($(FORTPATH),$(EMPTY))
  FORT=gfortran
  FORTPATH = $(strip $(shell which $(FORT) 2>/dev/null))
endif
ifeq ($(FORTPATH),$(EMPTY))
  FORT=sunf95
  FORTPATH = $(strip $(shell which $(FORT) 2>/dev/null))
endif
ifeq ($(FORTPATH),$(EMPTY))
  FORT=pgf95
  FORTPATH = $(strip $(shell which $(FORT) 2>/dev/null))
endif


FORT_NOPATH = $(notdir $(FORT))


######################################################################
##### FORTRAN COMPILER: Intel (ifort)                               ##
######################################################################
ifeq ($(FORT_NOPATH),ifort)
  FORTVERS=`$(FORT) -V 2>&1|head -1`
#  FLAGSF90=-traceback -check noarg_temp_created -check all
  FLAGSF90=-traceback -check noarg_temp_created
  FLAGSF77=-traceback
#  FLAGSFORALL=-align commons -warn noalignments -cpp -noD -vec-report-0 -msse2
  FLAGSFORALL=-align commons -warn noalignments -cpp -noD -msse2

# for version 12.1:
  FLAGSFORALL += -diag-disable 8290 -diag-disable 8291

  ifeq ($(STATIC),1)
    FLAGSFORALL += -static
  endif
  ifeq ($(FPE),3)
    FLAGSFORALL += -fpe3
  endif
  ifeq ($(FPE),2)
    FLAGSFORALL += -fpe2
  endif
  ifeq ($(FPE),1)
    FLAGSFORALL += -fpe1
  endif
  ifeq ($(FPE),0)
    FLAGSFORALL += -fpe0
  endif
  FLAGSDOUBLE=-r8
endif
######################################################################
##### FORTRAN COMPILER: GCC/gfortran (requires version 4.6+)        ##
######################################################################
ifeq ($(findstring gfortran,$(FORT_NOPATH)),gfortran)
  FORTVERS=`$(FORT) -v 2>&1 | grep -i 'gcc.version'`

  FLAGSF77 =-g -fbacktrace -fcheck=all
  FLAGSF77 += -Wall
  FLAGSF77 += -Wno-unused-variable
  FLAGSF77 += -Wno-unused-dummy-argument
  FLAGSF77 += -Wno-conversion
  FLAGSF77 += -Wno-maybe-uninitialized
  FLAGSF77 += -Wno-unused-label

  FLAGSF90 =-g -fbacktrace -fcheck=all
  FLAGSF90 += -Wall
  FLAGSF90 += -Wextra
#  FLAGSF90 += -Wuninitialized

  FLAGSFORALL=-ffree-line-length-none -Wno-align-commons
  ifeq ($(STATIC),1)
    FLAGSFORALL += -static -static-libgcc
  endif
  FLAGSDOUBLE=-fdefault-real-8 -fdefault-double-8
  ifeq ($(FPE),0)
#    FLAGSFORALL += -ffpe-trap=invalid,zero,overflow,underflow
    FLAGSFORALL += -ffpe-trap=invalid,zero,overflow,underflow,denormal
  endif
  ifeq ($(ARCH),32)
    FLAGSFORALL += -m32
  endif
endif
######################################################################
##### FORTRAN COMPILER: Sun/Oracle (sunf95, sunf90, ...)            ##
######################################################################
ifeq ($(findstring sunf,$(FORT_NOPATH)),sunf)
  FORTVERS=`$(FORT) -V 2>&1|head -1`
  FLAGSF90=-g -w2 -xcheck=%all
  FLAGSF77=$(FLAGSF90)
  ifeq ($(STATIC),1)
    FLAGSFORALL=-Bstatic
  endif
  FLAGSDOUBLE=-xtypemap=real:64
  ifeq ($(FPE),3)
    FLAGSFORALL += -ftrap=%none
  endif
  ifeq ($(FPE),1)
    FLAGSFORALL += -ftrap=common
  endif
  ifeq ($(FPE),0)
    FLAGSFORALL += -ftrap=%all
  endif
endif
######################################################################
##### FORTRAN COMPILER: Portland/PGI (pgfortran, pgf95, pgf90, ...) ##
######################################################################
ifeq ($(findstring pgf,$(FORT_NOPATH)),pgf)
  FORTVERS=`$(FORT) -V 2>&1|head -2`
  FLAGSF90=-g
  FLAGSF77=$(FLAGSF90)
  FLAGSDOUBLE=-r8
  FLAGSFORALL=-Mextend
  ifeq ($(STATIC),1)
    FLAGSFORALL += -Bstatic
  endif
endif
######################################################################
##### FORTRAN COMPILER: LLVM/flang                                  ##
######################################################################
ifeq ($(findstring flang,$(FORT_NOPATH)),flang)
  FORTVERS=`$(FORT) -v 2>&1|head -2`
  FLAGSF90=-g
  FLAGSF77=$(FLAGSF90)
  FLAGSDOUBLE=-r8
  FLAGSFORALL=-Mextend
  ifeq ($(STATIC),1)
    FLAGSFORALL += -Bstatic -static-flang-libs
  endif
endif

#########################################################
# Note: The compilers below are currently not supported.
#########################################################

#########################################################
##### FORTRAN COMPILER: g95
#########################################################
ifeq ($(FORT_NOPATH),g95)
  FORTVERS=`$(FORT) --version|head -1`
  FLAGSF90=-Wall -fbounds-check -ftrace=full
  FLAGSF77=$(FLAGSF90)
  FLAGSFORALL=-cpp -ffree-line-length-huge
  ifeq ($(STATIC),1)
    FLAGSFORALL += -static
  endif
  FLAGSDOUBLE=-r8
endif
#########################################################
##### FORTRAN COMPILER: Open64
#########################################################
ifeq ($(FORT_NOPATH),openf95)
  FORTVERS=`$(FORT) -V 2>&1|head -1`
  FLAGSF90=-g
  FLAGSF77=$(FLAGSF90)
  FLAGSDOUBLE=-r8
  ifeq ($(STATIC),1)
    FLAGSFORALL=--static -static-libgcc
  endif
endif
#########################################################
##### FORTRAN COMPILER: PathScale pathf95
#########################################################
ifeq ($(FORT_NOPATH),pathf95)
  FORTVERS=`$(FORT) -v 2>&1|head -1`
  FLAGSF90=-g -Wall -ffortran-bounds-check
  FLAGSF77=$(FLAGSF90)
  FLAGSDOUBLE=-r8
  ifeq ($(STATIC),1)
    FLAGSFORALL=-static
  endif
endif
#########################################################
##### FORTRAN COMPILER: Lahey lfc
#########################################################
ifeq ($(findstring lfc,$(FORT_NOPATH)),lfc)
  FORTVERS=`$(FORT) lfc --version  2>&1|head -1`
  # -chkglobal would be nice, however --chk x does not work with Pythia
  # FLAGSF90= --chk --chkglobal
  FLAGSF90= --chk
  FLAGSF77=
  FLAGSFORALL=--static --staticlink -Cpp
  FLAGSDOUBLE=-CcdRR8
endif
#########################################################
##### FORTRAN COMPILER: NAGWare nagfor
#########################################################
ifeq ($(FORT_NOPATH),nagfor)
  FORTVERS=`$(FORT) -v 2>&1|head -1`
  FLAGSF90=-g
  FLAGSF77=$(FLAGSF90)
#  FLAGSDOUBLE=-r8
  FLAGSDOUBLE=
  FLAGSFORALL=-Bstatic -f2003 -kind=byte -dusty -132
endif
#########################################################
##### FORTRAN COMPILER: Absoft f95
#########################################################
ifeq ($(FORT_NOPATH),f95)
  FORTVERS=`$(FORT) -v | head -1`
  FLAGSF90=
  FLAGSF77=
  FLAGSDOUBLE=-N113
  FLAGSFORALL=-TENV:simd_imask=on,simd_dmask=on,simd_zmask=on,simd_omask=on,simd_umask=on,simd_pmask=on -Rb -Rs -Rp  -W132 -m32 -p$(OBJDIR)
endif

#########################################################
##### OPTIMIZATION FLAGS (same for all compilers)
#########################################################
ifeq ($(MODE),opt)
  FLAGSF90=-O3
  FLAGSF77=-O3
endif
ifeq ($(MODE),opt0)
  FLAGSF90+=-O0
  FLAGSF77+=-O0
endif
ifeq ($(MODE),opt1)
  FLAGSF90=-O1
  FLAGSF77=-O1
endif
ifeq ($(MODE),opt2)
  FLAGSF90=-O2
  FLAGSF77=-O2
endif
ifeq ($(MODE),opt3)
  FLAGSF90=-O3
  FLAGSF77=-O3
endif
ifeq ($(MODE),opt4)
  # only works with gfortran & g95
  FLAGSF90=-O3 -ffast-math -funroll-loops -ftree-vectorize -march=native
  FLAGSF77=-O3 -ffast-math -funroll-loops -ftree-vectorize -march=native
endif
ifeq ($(MODE),opt5)
  # auto-parallelization for ifort & gfortran (with 4 threads)
  ifeq ($(FORT_NOPATH),ifort)
    FLAGSF90=-O3 -parallel -par-report1 -par-num-threads=4
    FLAGSF77=-O3 -parallel -par-report1 -par-num-threads=4
  endif
  ifeq ($(findstring gfortran,$(FORT_NOPATH)),gfortran)
    FLAGSF90=-O3 -ftree-parallelize-loops=4
    FLAGSF77=-O3 -ftree-parallelize-loops=4
  endif
endif

# for gfortran:
ifeq ($(MODE),lto)
  FLAGSF90=-O3 -flto=3 -fuse-linker-plugin
  FLAGSF77=-O3 -flto=3 -fuse-linker-plugin
  AR = gcc-ar
endif

#for ifort:
ifeq ($(MODE),ipo)
# following 2 lines are connected! if -ipo, then also 'xiar' instead of 'ar'
  FLAGSFORALL += -O3 -ipo
  AR = xiar
endif

#########################################################
##### PROFILING FLAGS (for ifort & gfortran)
#########################################################
ifeq ($(MODE),prof)
  FLAGSF90+= -pg
  FLAGSF77+= -pg
endif
#########################################################
##### FLAGS for producing a Call Graph (gfortran only)
#########################################################
ifeq ($(MODE),callGraph)
  FLAGSFORALL+=-fdump-rtl-expand
endif
#########################################################

### Adding some flags given at the command line
FLAGSFORALL += $(ARGS)


export FORT
export FORT_NOPATH
export FORTVERS
export FLAGSF90
export FLAGSF77
export FLAGSFORALL
export FLAGSDOUBLE

#*******************************************************************************
#*******************************************************************************
#************* 3) Targets
#*******************************************************************************
#*******************************************************************************


#*******************************************************************************
#****e* Makefile/all
# NAME
# make all
# PURPOSE
# This is the main target. It builds the GiBUU executable.
#*******************************************************************************
.PHONY : all
all: code/Makefile
	@$(MAKE) $(noPrintDirectory) doCompile
	@$(MAKE) $(noPrintDirectory) buildExe


# automatically make renew
code/Makefile:
	@$(MAKE) $(noPrintDirectory) renew


#*******************************************************************************
#****e* Makefile/bw
# NAME
# make bw
# PURPOSE
# Black and white compiling. Useful for piping compilation output to a log file.
#*******************************************************************************
.PHONY : bw
bw:
	@$(MAKE) StartHeader= EndHeader= StartHeader_green= StartHeader_red= blue= endBlue=


.PHONY : quick
quick:
	@$(MAKE) cleanEXE
	@$(MAKE) doCompile
	@$(MAKE) buildExe

.PHONY : printCompiler
printCompiler:
	@echo "Operating System:" $(OS_LONG)
ifeq ($(FORTPATH),$(EMPTY))
	@echo "ERROR: No supported compiler found!"
else
	@echo "The compiler " $(FORT_NOPATH) " is used (" $(FORTPATH) ")."
	@echo $(FORTVERS)
endif

.PHONY : print
print:
	@echo "!!! CURDIR     =" $(CURDIR)
	@echo "!!! SUBDIR     =" $(SUBDIR)
	@echo "!!! SRC        =" $(SRC)

.PHONY : subdirs $(SUBDIR)
subdirs:
	@$(ECHO) $(StartHeader)Collecting all source code information...$(EndHeader)
	@for X in $(SUBDIR); do \
	  (cd $$X && $(MAKE) $(noPrintDirectory) iterate;)\
	done
	@$(MAKE) printFinished $(noPrintDirectory)

.PHONY : printFinished
printFinished:
	@$(ECHO) $(StartHeader)...finished$(EndHeader)

.PHONY : doCompile
doCompile:
	@$(MAKE) printCompiler $(noPrintDirectory)
	@$(MAKE) subdirs $(noPrintDirectory)
	@scripts/removeDeadLinks.sh
	@$(ECHO) $(StartHeader)Logging code version...$(EndHeader)
	@$(MAKE) $(noPrintDirectory) version
	@$(MAKE) printFinished $(noPrintDirectory)
	@$(ECHO) $(StartHeader)Compiling source code...$(EndHeader)
	@cd objects && $(MAKE) $(noPrintDirectory) compileOBJ
	@$(MAKE) printFinished $(noPrintDirectory)
	@$(ECHO) $(StartHeader)Building libraries...$(EndHeader)
	@cd objects && $(MAKE) $(noPrintDirectory) buildLIBS
	@$(MAKE) printFinished $(noPrintDirectory)

.PHONY : doCompile
buildExe:
	@cd objects && $(MAKE) $(noPrintDirectory) EXE
	@ln -sf $(OBJDIR)/$(NameOfExe).x $(EXEDIR)/$(NameOfExe).x
	@if [ -f $(EXEDIR)/$(NameOfExe).x ]; then \
	  $(ECHO) $(StartHeader_green)SUCCESS: $(NameOfExe).x generated.$(EndHeader) ; \
	else \
	  $(ECHO) $(StartHeader_red)!!!!! ERROR: $(NameOfExe).x not generated !!!!!$(EndHeader) ; \
	fi


#*******************************************************************************
#****e* Makefile/clean
# NAME
# make clean
# PURPOSE
# Removes all "#*" and "*~" files from the tree. These are files generated
# by editing the code.
#*******************************************************************************
.PHONY : clean
clean:
	@$(ECHO) $(StartHeader)Cleaning up ...$(EndHeader)
	@$(FIND) ./ -name "#*" -exec rm '{}' ';'
	@$(FIND) ./ -name "*~" -exec rm '{}' ';'


#*******************************************************************************
#****e* Makefile/cleanEXE
# NAME
# make cleanEXE
# PURPOSE
# Remove all "*.x" files from the tree.
#*******************************************************************************
.PHONY : cleanEXE
cleanEXE:
	@$(FIND) ./code -name "*.x" -exec rm '{}' ';'
	@$(FIND) ./testRun -name "*.x" -exec rm '{}' ';'
	@$(FIND) ./objects -name "*.x" -exec rm '{}' ';'


#*******************************************************************************
#****e* Makefile/cleanOBJ
# NAME
# make cleanOBJ
# PURPOSE
# Remove all "*.o" and "*.mod" files from the directory "objects".
#*******************************************************************************
.PHONY : cleanOBJ
cleanOBJ:
	@-cd objects && rm -f *.o *.mod

#*******************************************************************************
#****e* Makefile/cleanLIB
# NAME
# make cleanLIB
# PURPOSE
# Remove all "*.a"  files from the directory "objects/LIB/lib".
#*******************************************************************************
.PHONY : cleanLIB
cleanLIB:
	@-cd objects/LIB/lib && rm -f *.a


#*******************************************************************************
#****e* Makefile/veryclean
# NAME
# make veryclean
# PURPOSE
# First calls the targets "clean", "cleanEXE" and "cleanOBJ",
# then removes all "*.o" and "*.mod" files from the tree, removes all symlinks
# in the objects/ directory and deletes version.f90.
#*******************************************************************************
.PHONY : veryclean
veryclean: clean cleanEXE cleanOBJ
	@$(FIND) ./ -name "*.o" -exec rm '{}' ';'
	@$(FIND) ./ -name "*.mod" -exec rm '{}' ';'
	@rm -f objects/*.f90
	@rm -f objects/*.f
	@rm -f objects/*.F90
	@rm -f objects/*.F
	@rm -f objects/*.expand
	@rm -f objects/*.inc
	@rm -f code/inputOutput/version.f90


#*******************************************************************************
#****e* Makefile/superclean
# NAME
# make superclean
# PURPOSE
# First calls the target "veryclean",
# then removes all "fort.*" files from the tree.
#*******************************************************************************
.PHONY : superclean
superclean: veryclean
	@$(FIND) ./ -name "fort.????" -exec rm '{}' ';'
	@$(FIND) ./ -name "fort.???" -exec rm '{}' ';'
	@$(FIND) ./ -name "fort.??" -exec rm '{}' ';'
	@$(FIND) ./ -name "fort.?" -exec rm '{}' ';'


#*******************************************************************************
#****e* Makefile/svnclean
# NAME
# make svnclean
# PURPOSE
# Remove all files which are not included in the svn repository.
# NOTES
# This is very dangerous, since it removes all unversioned files in the
# GiBUU directory! Don't use it unless you know exactly what you're doing!
#*******************************************************************************
.PHONY : svnclean
svnclean:
	@svn status | grep ^? | awk '{print $2}' | xargs rm -rf


#*******************************************************************************
#****e* Makefile/jobcard
# NAME
# make jobcard
# PURPOSE
# This target builds a master job card in testRun/jobCards,
# which includes all switches.
#*******************************************************************************
.PHONY : jobcard
jobcard:
	@cd testRun/jobCards && . makeMasterJobCard


#*******************************************************************************
#****e* Makefile/doku
# NAME
# make doku
# PURPOSE
# This target rebuilds the whole html documentation by calling ROBODoc.
#*******************************************************************************
.PHONY : doku
doku:
	robodoc --rc robodoc/GenerateDoku.rc
	robodoc --rc robodoc/MakefileDoku.rc
	@cp robodoc/index.html.head ../Documentation/index.html
	@cp robodoc/robodoc.css     ../Documentation/robodoc.css
	@echo "  Generated from release 2017 `svnversion -n .` with ROBODoc `robodoc --version` on `date` "  >> ../Documentation/index.html
	@echo "</div>"  >> ../Documentation/index.html
	@echo "</body>" >> ../Documentation/index.html
	@echo "</html>" >> ../Documentation/index.html


#*******************************************************************************
#****e* Makefile/showDoku
# NAME
# make showdoku
# PURPOSE
# This target opens the documentation with firefox.
#*******************************************************************************
.PHONY : showDoku
showDoku:
	@echo "Opening the HTML documentation..."
	@firefox ../Documentation/index.html
	@echo "... done!"


#*******************************************************************************
#****e* Makefile/svnupdate
# NAME
# make svnupdate
# PURPOSE
# This target performs an svn update of the code and input.
#*******************************************************************************
.PHONY : svnupdate
svnupdate:
	@$(ECHO) 'Updating source code via SVN'
	@svn update
	@$(ECHO) 'Updating input files via SVN'
	@if [ -d $(PATH_TO_INPUT) ]; then \
	  svn update $(PATH_TO_INPUT) ;\
	else \
	  $(ECHO) "WARNING: The buuinput directory is not available in its usual path.";\
	  $(ECHO) "         => Updating not possible!";\
	fi


#*******************************************************************************
#****e* Makefile/update
# NAME
# make update
# PURPOSE
# This target performs an svn update and thereafter recompiles the whole code.
#*******************************************************************************
.PHONY : update
update: svnupdate
	$(MAKE) renew
	$(MAKE) all


#*******************************************************************************
#****e* Makefile/allTests
# NAME
# make allTests
# PURPOSE
# This target builds all test cases
# (the testRun directories are listed in "Makefile.testRun.List").
#*******************************************************************************
.PHONY : allTests
allTests:
	@$(MAKE) doCompile
	@for X in `cat Makefile.testRun.List|grep -v "#"`; do \
	  cd $(CURDIR) && cd $$X && $(MAKE) iterate; \
	done

#*******************************************************************************

.PHONY : CalledFromTestRunFull CalledFromTestRun

CalledFromTestRunFull:
#@echo $(CALLINGDIR)
	@$(MAKE) doCompile
	@cd $(CALLINGDIR) && $(MAKE) iterate;

CalledFromTestRun:
#@echo $(CALLINGDIR)
	@cd $(CALLINGDIR) && $(MAKE) iterate;

#*******************************************************************************

### ONLY for MAKEFILE-Maintainers:
.PHONY : Makefiles MAKEFILEclean
Makefiles: #print
	@$(ECHO) $(StartHeader)Setting up Makefiles ...$(EndHeader)
	@for X in $(SUBDIR); do\
	  (cp $(ROOTDIR)/Makefile.SUBlink $$X/Makefile; \
	  cd $$X && $(MAKE) Makefiles;)\
	done
	@sed -e s\|WORKDIR\|"$(ROOTDIR)"\| $(ROOTDIR)/Makefile.objects > mmm.2
	@cp mmm.2 objects/Makefile
	@touch objects/MakefileDepend
	@for X in `cat Makefile.testRun.List|grep -v "#"`; do \
	  cp mmm.2 $$X/Makefile; \
	  touch $$X/MakefileDepend; \
	done
	@rm -f mmm.2

MAKEFILEclean: veryclean
#	$(FIND) ./ -name "Makefile" -exec rm '{}' ';'
	@$(FIND) ./ -name "MakefileDepend" -exec rm '{}' ';'


#*******************************************************************************
#****e* Makefile/renew
# NAME
# make renew
# PURPOSE
# If dramatical changes were done on the tree structure,
# directories added or deleted etc you should call
# "make renew" before you call "make".
#
# This target implies a call of "veryclean". Then
# in every directory the corresponding local Makefile is
# deleted and finally updated.
#*******************************************************************************

.PHONY : renew

#renew : input MAKEFILEclean Makefiles
renew : MAKEFILEclean Makefiles


#*******************************************************************************
#****e* Makefile/cleanTestRun
# NAME
# make cleanTestRun
# PURPOSE
# Clean up the directory "testRun" (by removing files like *.dat, *.tex, fort.*,
# main.run and waitBar.eps).
#*******************************************************************************
.PHONY : cleanTestRun
cleanTestRun:
	-rm -rf $(ROOTDIR)/testRun/*.dat
	-rm -rf $(ROOTDIR)/testRun/*.tex
	-rm -rf $(ROOTDIR)/testRun/fort.*
	-rm -rf $(ROOTDIR)/testRun/main.run
	-rm -rf $(ROOTDIR)/testRun/waitBar.eps



#*******************************************************************************
#****e* Makefile/CallGraph
# NAME
# make CallGraph
# PURPOSE
# Create call-graph stuff.
# NOTES
# This target is not working from scratch.
# You have to call "make renew" and then "make CallGraph FORT=gfortran"
# or similar, since only gfortran reports the necessary RTL files.
#*******************************************************************************
.PHONY : CallGraph
CallGraph:
	make MODE=callGraph
#	egypt objects/*.expand > egypt.out
	scripts/Own_Egypt_CallGraph.pl objects/*.expand > Own_Egypt_CallGraph.out


#*******************************************************************************
#****e* Makefile/PlotCallGraph
# NAME
# make PlotCallGraph
# PURPOSE
# Creates a plot of the call graph in svg format.
# NOTES
# Requires the 'graphviz' tools (www.graphviz.org/). The perl script 'egypt'
# v1.10 is used (www.gson.org/egypt/). Only works with gfortran.
#*******************************************************************************
.PHONY : PlotCallGraph
PlotCallGraph:
	make MODE=callGraph
	time scripts/egypt.pl objects/*.expand > egypt.out
	time dot -Grankdir=LR -Tsvg -o CallGraph.svg egypt.out


#*******************************************************************************
#****e* Makefile/ModGraph
# NAME
# make ModGraph
# PURPOSE
# Creates a plot of the module dependency graph in svg format.
# NOTES
# Requires the 'graphviz' tools (www.graphviz.org/). Only works with gfortran.
#*******************************************************************************
.PHONY : ModGraph
ModGraph:
	@$(MAKE) subdirs $(noPrintDirectory)
	@scripts/removeDeadLinks.sh
	@$(MAKE) $(noPrintDirectory) version
	@cd objects && $(MAKE) $(noPrintDirectory) ModGraph


#*******************************************************************************
#****e* Makefile/version
# NAME
# make version
# PURPOSE
# Renew version info.
#*******************************************************************************
.PHONY : version
version: objects/version.f90
#	@$(MAKE) objects/version.f90

objects/version.f90 : $(allSrcFiles)
	@echo 'module version'                                           >  ./code/inputOutput/version.f90
	@echo ' contains  '                                              >> ./code/inputOutput/version.f90
	@echo ' subroutine PrintVersion '                                >> ./code/inputOutput/version.f90
	@echo ' implicit none '                                          >> ./code/inputOutput/version.f90
	@echo ' write(*,2) "Release     : 2017" '                        >> ./code/inputOutput/version.f90
	@echo ' write(*,2) "SVN Revision: '`svnversion -n .` '" '        >> ./code/inputOutput/version.f90
	@echo ' write(*,2) "Build Date  : '`date` '" '                   >> ./code/inputOutput/version.f90
	@echo ' write(*,2) "OS          : '$(OS_LONG)'" '                >> ./code/inputOutput/version.f90
	@echo ' write(*,2) "Compiler    : '$(FORTVERS)'" '               >> ./code/inputOutput/version.f90
	@echo ' write(*,2) "Flags       : $(FLAGSF90) $(FLAGSFORALL)" '  >> ./code/inputOutput/version.f90
	@echo ' write(*,2) "PATH        : '`pwd -P` '" '                 >> ./code/inputOutput/version.f90
	@svn status -q |awk -F @ '{printf " write(*,2) \" %s\"\n", $$1}' >> ./code/inputOutput/version.f90
	@echo ' 2 FORMAT(" ",A)'                                         >> ./code/inputOutput/version.f90
	@echo ' end subroutine PrintVersion '                            >> ./code/inputOutput/version.f90
	@echo 'end module version '                                      >> ./code/inputOutput/version.f90
	@ln -sf $(ROOTDIR)/code/inputOutput/version.f90 $(ROOTDIR)/objects/version.f90


#####################################################################

.PHONY : buildLHAPDF buildPDFLIB buildPDFstubb buildPDF

buildLHAPDF: printCompiler
	@cd $(LIBDIRSRC)/PDF && $(MAKE) LHAPDF;

buildPDFLIB: printCompiler
	@cd $(LIBDIRSRC)/PDF && $(MAKE) PDFLIB;

buildPDFstubb: printCompiler

buildPDF:

.PHONY : buildBZIP2_32 buildBZIP2

buildBZIP2_32: printCompiler
	@cd $(LIBDIRSRC)/BZ2 && $(MAKE) LIBBZ2_32;

buildBZIP2: printCompiler
	@cd $(LIBDIRSRC)/BZ2 && $(MAKE) LIBBZ2;

# DO NOT DELETE
