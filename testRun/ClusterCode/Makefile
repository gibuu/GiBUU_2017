######################################################### -*- Makefile -*-
#****e* /Makefile
# NAME
# Makefile
# PURPOSE
# This is the main "Makefile" of the ClusterCode. 
# NOTES
# Copied from the GiBUU workingCode.
#***
#########################################################


#********************************************************
#********************************************************
#************* 1) Setting variables
#********************************************************
#********************************************************

.SUFFIXES:            # Delete the default suffixes

### PROGRAMS:

export SHELL = /bin/bash
export OS_LONG = '$(shell uname -s -r -m)'
export COPY=-cp -p -u
export noPrintDirectory=--no-print-directory
export ECHO = echo -e
export PERL = $(strip $(shell which perl 2>/dev/null))
ifeq ($(PERL),$(EMPTY))
  export MAKEDEP = $(strip $(shell which makedepf90 2>/dev/null))
endif

### DIRECTORIES:

export ROOTDIR := $(CURDIR)

# if CURDIR is a symbolc link, then use better the following:
export ROOTDIRLONG := `pwd -P`

# Directory to store the *.o and *.mod files:
export OBJDIR := $(CURDIR)/objects

# Directory to store the executable:
export EXEDIR := $(CURDIR)/ClusterRuns

## Directory to find sources of additional libraries:
#export LIBDIRSRC := $(ROOTDIRLONG)/../libraries

## Directory to find the .a files:
#export LIBDIR := $(OBJDIR)/LIB/lib


# The variable "allSrcFiles" includes :
# * path of all the source files which are linked to objects
# * path to the main.f90 file. 
SRCf90 := $(wildcard objects/*.f90)
SRCF90 := $(wildcard objects/*.F90)
SRCf77 := $(wildcard objects/*.f)
allSrcFiles:= $(SRCf90) $(SRCF90) $(SRCf77) code/main.f90
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
#********************************************************
#********************************************************
#************* 2) Defining the compiler
#********************************************************
#********************************************************

EMPTY =
FPE = 3
STATIC = 0

### default: Intel

FORT=ifort
FORTPATH = $(strip $(shell which $(FORT) 2>/dev/null))

### fallback options if ifort is not available

ifeq ($(FORTPATH),$(EMPTY))
  FORT=gfortran
  FORTPATH = $(strip $(shell which $(FORT) 2>/dev/null))
endif

ifeq ($(FORTPATH),$(EMPTY))
  FORT=gfortran-4
  FORTPATH = $(strip $(shell which $(FORT) 2>/dev/null))
endif

ifeq ($(FORTPATH),$(EMPTY))
  FORT=g95
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


#########################################################
##### FORTRAN COMPILER: INTEL
#########################################################
ifeq ($(FORT_NOPATH),ifort)
  FORTVERS=`$(FORT) -V 2>&1|head -1`
  FLAGSF90=-traceback -check noarg_temp_created -check all
  FLAGSF77=-traceback
  FLAGSFORALL=-align commons -warn noalignments -cpp -noD -vec-report-0 -msse2

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
#########################################################
##### FORTRAN COMPILER: GCC/gfortran
##### (requires GCC 4.3.0+)
#########################################################
ifeq ($(findstring gfortran,$(FORT_NOPATH)),gfortran)
  FORTVERS=`$(FORT) -v 2>&1|tail -1`
  FLAGSF90=-Wall -Wtabs -fbounds-check -g -fbacktrace
  FLAGSF77=$(FLAGSF90)
  FLAGSFORALL=-ffree-line-length-none
  ifeq ($(STATIC),1)
    FLAGSFORALL += -static -static-libgcc
  endif
  FLAGSDOUBLE=-fdefault-real-8 -fdefault-double-8
  ifeq ($(FORT_NOPATH),gfortran-4.4)
    FLAGSFORALL += -Wno-align-commons
  endif
  ifeq ($(FORT_NOPATH),gfortran-4.5)
    FLAGSFORALL += -Wno-align-commons
    FLAGSF90 += -fcheck=all
    FLAGSF77=$(FLAGSF90)
  endif
  ifeq ($(FPE),0)
    FLAGSFORALL += -ffpe-trap=invalid,zero,overflow,underflow
  endif
  ifeq ($(ARCH),32)
    FLAGSFORALL += -m32
  endif
endif
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
##### FORTRAN COMPILER: Sun/Oracle
#########################################################
ifeq ($(findstring sunf95,$(FORT_NOPATH)),sunf95)
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
#########################################################
##### FORTRAN COMPILER: Portland pgf95
#########################################################
ifeq ($(FORT_NOPATH),pgf95)
  FORTVERS=`$(FORT) -V 2>&1|head -2`
  FLAGSF90=-g
  FLAGSF77=$(FLAGSF90)
  FLAGSDOUBLE=-r8
  FLAGSFORALL=-Mextend
  ifeq ($(STATIC),1)
    FLAGSFORALL += -Bstatic
  endif
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
  FLAGSDOUBLE=-r8
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
ifeq ($(MODE),opt4)   # only works with gfortran & g95
  FLAGSF90=-O3 -ffast-math -funroll-loops -ftree-vectorize -march=native
  FLAGSF77=-O3 -ffast-math -funroll-loops -ftree-vectorize -march=native
endif
#########################################################
##### PROFILING FLAGS (for ifort & gfortran)
#########################################################
ifeq ($(MODE),prof)
  FLAGSF90+= -pg
  FLAGSF77+= -pg
endif
#########################################################

export FORT
export FORT_NOPATH
export FORTVERS
export FLAGSF90
export FLAGSF77
export FLAGSFORALL
export FLAGSDOUBLE

#********************************************************
#********************************************************
#************* 3) Targets
#********************************************************
#********************************************************

#****e* Makefile/all
# NAME
# make all
# PURPOSE
# This is the first and main target: this builds the executable "main.x"
#***
all: code/Makefile
	@$(MAKE) $(noPrintDirectory) doCompile
	@$(MAKE) $(noPrintDirectory) buildExe


# automatically make renew
code/Makefile:
	@$(MAKE) $(noPrintDirectory) renew


# black and white compiling
bw: 
	@$(MAKE) StartHeader= EndHeader= StartHeader_green= StartHeader_red= blue= endBlue= 


quick:
	@$(MAKE) cleanEXE
	@$(MAKE) doCompile
	@$(MAKE) buildExe


printCompiler:
	@echo "Operating System:" $(OS_LONG)
ifeq ($(FORTPATH),$(EMPTY))
	@echo "ERROR: No supported compiler found!"
else
	@echo "The compiler " $(FORT_NOPATH) " is used (" $(FORTPATH) ")."
	@echo $(FORTVERS)
endif 

#iterate: printCompiler subdirs doCompile
#	@if [ -f $(EXEDIR)/main.x ]; then \
#	$(ECHO) $(StartHeader_green)SUCCESS: main.x generated.$(EndHeader) ; \
#	else \
#	$(ECHO) $(StartHeader_red)!!!!! ERROR: main.x not generated !!!!!$(EndHeader) ; \
#	fi

print:
	@echo "!!! CURDIR     =" $(CURDIR)
	@echo "!!! SUBDIR     =" $(SUBDIR)
	@echo "!!! SRC        =" $(SRC)

subdirs:
	@$(ECHO) $(StartHeader)Collecting all source code information...$(EndHeader)
	@for X in $(SUBDIR); do \
	  (cd $$X && $(MAKE) $(noPrintDirectory) iterate;)\
	done	
	@$(MAKE) printFinished $(noPrintDirectory)


printFinished:
	@$(ECHO) $(StartHeader)...finished$(EndHeader)


doCompile:
	@$(MAKE) printCompiler $(noPrintDirectory)
	@$(MAKE) subdirs $(noPrintDirectory)
	@../../scripts/removeDeadLinks.sh
	@$(ECHO) $(StartHeader)Logging code version...$(EndHeader)
	@$(MAKE) $(noPrintDirectory) version
	@$(MAKE) printFinished $(noPrintDirectory)
	@$(ECHO) $(StartHeader)Compiling source code...$(EndHeader)
	@cd objects && $(MAKE) $(noPrintDirectory) compileOBJ
	@$(MAKE) printFinished $(noPrintDirectory)
#	@$(ECHO) $(StartHeader)Building libraries...$(EndHeader)
#	@cd objects && $(MAKE) $(noPrintDirectory) buildLIBS
#	@$(MAKE) printFinished $(noPrintDirectory)
#	@$(ECHO) $(StartHeader)Building executable...$(EndHeader)

#	@cd objects && $(MAKE) $(noPrintDirectory) EXE
#	@ln -sf $(CURDIR)/objects/main.x $(CURDIR)/ClusterRuns/main.x
#	@$(MAKE) printFinished $(noPrintDirectory)


buildExe:
	@cd objects && $(MAKE) $(noPrintDirectory) EXE
	@ln -sf $(CURDIR)/objects/main.x $(CURDIR)/ClusterRuns/main.x
	@if [ -f $(CURDIR)/objects/main.x ]; then \
	  $(ECHO) $(StartHeader_green)SUCCESS: main.x generated.$(EndHeader) ; \
	else \
	  $(ECHO) $(StartHeader_red)!!!!! ERROR: main.x not generated !!!!!$(EndHeader) ; \
	fi


#****e* Makefile/cleanEXE
# NAME
# make cleanEXE
# PURPOSE
# remove every "*.x" files from the tree 
#***
cleanEXE:
	@find ./code -name "*.x" -exec rm '{}' ';'
	@find ./ClusterRuns -name "*.x" -exec rm '{}' ';'
	@find ./objects -name "*.x" -exec rm '{}' ';'

#****e* Makefile/cleanOBJ
# NAME
# make cleanOBJ
# PURPOSE
# remove every "*.o" and "*.mod" files from the 
# directory "objects"
#***
cleanOBJ:
	@-cd objects && rm -f *.o *.mod

#****e* Makefile/veryclean
# NAME
# make veryclean
# PURPOSE
# first calls the targets "cleanEXE" and "cleanOBJ",
# then removes every "*.o" and "*.mod" files from the tree.
#***
veryclean: cleanEXE cleanOBJ
	@find ./ -name "*~" -exec rm '{}' ';'
	@find ./ -name "*.o" -exec rm '{}' ';'
	@find ./ -name "*.mod" -exec rm '{}' ';'
	@rm -f objects/version.f90

#****e* Makefile/superclean
# NAME
# make superclean
# PURPOSE
# first calls the targets "veryclean" , 
# then removes all "#*", "fort.*" and all "*.dat" (besides the
# ones in buuinput) files from the tree.
#***
superclean: veryclean
	@find ./ -name "#*" -exec rm '{}' ';'
	@find ./ -name "fort.????" -exec rm '{}' ';'
	@find ./ -name "fort.???" -exec rm '{}' ';'
	@find ./ -name "fort.??" -exec rm '{}' ';'
	@find ./ -name "fort.?" -exec rm '{}' ';'



#****e* Makefile/update
# NAME
# make update
# PURPOSE
# This target performs svn update and thereafter recompiles the whole code.
#***
svnupdate:
	@echo 'Updating the code via SVN'
	@svn update

update: svnupdate
	$(MAKE) renew
	$(MAKE) all

CalledFromClusterRunsFull:
	@echo $(CALLINGDIR)
#	@echo $(ROOTDIR)
	$(MAKE)
	cd $(CALLINGDIR) && $(MAKE) iterate;

CalledFromClusterRuns:
	@echo $(CALLINGDIR)
#	@echo $(ROOTDIR)
	cd $(CALLINGDIR) && $(MAKE) iterate;


### ONLY for MAKEFILE-Maintainers:

Makefile.ClusterRuns.List:
	touch Makefile.ClusterRuns.List

MAKEFILE: Makefile.ClusterRuns.List #print
	@$(ECHO) $(StartHeader)Setting up Makefiles ...$(EndHeader)
	@for X in $(SUBDIR); do\
	  (cp $(ROOTDIR)/Makefile.SUBlink $$X/Makefile; \
	  cd $$X && $(MAKE) MAKEFILE;)\
	done
	@sed -e s\|WORKDIR\|"$(ROOTDIR)"\| $(ROOTDIR)/Makefile.objects > mmm.2
	@cp mmm.2 objects/Makefile
#	@touch objects/MakefileDepend
	@for X in `cat Makefile.ClusterRuns.List|grep -v "#"`; do \
	  cp mmm.2 $$X/Makefile; \
	done
	@rm -f mmm.2

MAKEFILEclean: veryclean
#	@find ./ -name "Makefile" -exec rm '{}' ';'
	@find ./ -name "MakefileDepend" -exec rm '{}' ';'


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
#***
renew : MAKEFILEclean MAKEFILE 


#****e* Makefile/cleanClusterRuns
# NAME
# make cleanClusterRuns
# PURPOSE
# clean up the directory "ClusterRuns"
#***
cleanClusterRuns:
	rm -rf $(ROOTDIR)/ClusterRuns/*.dat
	rm -rf $(ROOTDIR)/ClusterRuns/*.tex
	rm -rf $(ROOTDIR)/ClusterRuns/fort.*
	rm -rf $(ROOTDIR)/ClusterRuns/main.run
	rm -rf $(ROOTDIR)/ClusterRuns/waitBar.eps



#####################################################################
### Renew version info
#####################################################################
version:
	@$(MAKE) objects/version.f90

objects/version.f90 : $(allSrcFiles) 
	@echo 'module version'                                           >  ./code/version.f90 
	@echo ' contains  '                                              >> ./code/version.f90
	@echo ' subroutine PrintVersion '                                >> ./code/version.f90
	@echo ' implicit none '                                          >> ./code/version.f90
	@echo ' write(*,2) "SVN Revision: '`svnversion -n .` '" '        >> ./code/version.f90
	@echo ' write(*,2) "Build Date  : '`date` '" '                   >> ./code/version.f90
	@echo ' write(*,2) "OS          : '$(OS_LONG)'" '                >> ./code/version.f90
	@echo ' write(*,2) "Compiler    : '$(FORTVERS)'" '               >> ./code/version.f90
	@echo ' write(*,2) "Flags       : $(FLAGSF90) $(FLAGSFORALL)" '  >> ./code/version.f90
	@echo ' write(*,2) "PATH        : '`pwd -P` '" '                 >> ./code/version.f90
	@svn status -q |awk -F @ '{printf " write(*,2) \" %s\"\n", $$1}' >> ./code/version.f90
	@echo ' 2 FORMAT(" ",A)'                                         >> ./code/version.f90
	@echo ' end subroutine PrintVersion '                            >> ./code/version.f90
	@echo 'end module version '                                      >> ./code/version.f90
	@ln -sf $(ROOTDIR)/code/version.f90 $(ROOTDIR)/objects/version.f90


#####################################################################



#buildLHAPDF: printCompiler
#	@cd $(LIBDIRSRC)/PDF && $(MAKE) LHAPDF;
#
#buildPDFLIB: printCompiler
#	@cd $(LIBDIRSRC)/PDF && $(MAKE) PDFLIB;
#
#buildPDFstubb: printCompiler
#
#buildPDF:

# DO NOT DELETE
