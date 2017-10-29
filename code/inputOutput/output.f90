!******************************************************************************
!****m* /output
! NAME
! module output
! PURPOSE
! Includes basic tools for output.
!******************************************************************************
module output
  implicit none
  private


  !****************************************************************************
  !****g* output/line
  ! PURPOSE
  ! present a string of length 79 with '*'
  ! (not 80 because of f90-internal line breaking!!!)
  ! SOURCE
  character(79), parameter, public :: line='*****************************************************************'
  !****************************************************************************


  !****************************************************************************
  !****g* output/header
  ! PURPOSE
  ! Format of the BUU-Header in the output
  ! USAGE
  ! write(*,header) 'Hallo'
  ! SOURCE
  character(100), parameter, public :: header='(79("#"),/,79("#"),/,10("#")," ",A,/,79("#"),/,79("#"),/)'
  !****************************************************************************


  !****************************************************************************
  !****g* output/chapter
  ! PURPOSE
  ! Format of the output announcing the start of some chapter of the code.
  ! USAGE
  ! write(*,chapter) 'Hallo'
  ! SOURCE
  character(100), parameter, public :: chapter='(79("#"),/,5("#")," ",A,/,79("#"),/)'
  !****************************************************************************


  !****************************************************************************
  !****g* output/subchapter
  ! PURPOSE
  ! Format of the output announcing the start of some subchapter of the code.
  ! USAGE
  ! write(*,subchapter) 'Hallo'
  ! SOURCE
  character(100), parameter, public :: subchapter='(79("-"),/,5("-")," ",A,/,79("-"),/)'
  !****************************************************************************


  !****************************************************************************
  !****g* output/paragraph
  ! PURPOSE
  ! Format of the output announcing the start of some paragraph of the code.
  ! USAGE
  ! write(*,paragraph) 'Hallo'
  ! SOURCE
  character(100), parameter, public :: paragraph='(5("~")," ",A," ",5("~"))'
  !****************************************************************************


  integer :: iDoPrLevelDefault
  !****************************************************************************
  !****g* output/DoPrLevelDefault
  ! PURPOSE
  ! Default values for DoPrLevel
  ! SOURCE
  logical, dimension(-10:5), public, parameter :: DoPrLevelDefault = &
       (/ (.FALSE.,iDoPrLevelDefault=-10,-1),(.TRUE.,iDoPrLevelDefault=0,5)/)
  !****************************************************************************


  !****************************************************************************
  !****g* output/DoPrLevel
  ! PURPOSE
  ! Flags, whether message at (error)level n should be printed or not.
  ! The (error) levels are:
  ! * Level 5 (TERMINAL) -- A terminal error is not an informational error
  !   because corrective action within the program is generally not reasonable.
  !   In normal usage, execution should be terminated immediately when an
  !   error of this class occurs.
  ! * Level 4 (FATAL) -- A fatal error indicates the existence of a
  !   condition that may be serious. In most cases, user or calling
  !   routine must take corrective action to recover.
  ! * Level 3 (WARNING) -- A warning indicates the existence of a
  !   condition that may require corrective action by user or calling
  !   routine
  ! * Level 2 (ALERT) -- Indicates that the user should be adviced about
  !   events occuring in the code.
  ! * Level 1 (NOTE) -- Is issued to indicate the possibility of a trivial
  !   error or simply to provide information about the computations.
  ! * Level 0 (BASIC) -- Some basic messages
  ! The levels -10 up to -1 are left for the programmer for individual
  ! purposes. One can switch on/off writing of message independant of
  ! the error classification.
  ! They are used for...:
  ! *  -1: switching on/off the messages:
  !    "WARNING: PYEVNT-Loop -> failure",
  !    "DoColl_nuN: itry= ..."
  ! *  -2: ...(not used)
  ! * etc
  ! SOURCE
  logical, dimension(-10:5), public :: DoPrLevel = DoPrLevelDefault
  !****************************************************************************


  !****************************************************************************
  !****s* output/WriteParticle
  ! NAME
  ! subroutine WriteParticle(iFile,iEnsemble,iPart, Part)
  ! subroutine WriteParticle(iFile,iEnsemble, Part)
  !
  ! PURPOSE
  ! Write the stored information about a particle to a file/stdout.
  !
  ! INPUTS
  ! * integer        :: iFile     -- number of file (6=stdout)
  ! * integer        :: iEnsemble -- number of Ensemble (maybe arbitrary)
  ! * integer        :: iPart     -- number of Particle (maybe arbitrary)
  ! * type(particle) :: Part      -- Particle to print
  ! or:
  ! * integer        :: iFile     -- number of file (6=stdout)
  ! * integer        :: iEnsemble -- number of Ensemble
  ! * type(particle),dimension(:) :: Part -- Particles to print
  !
  ! NOTES
  ! * This is the format prefered by KG ;)
  ! * If the particles are given in an array, all non-zero entries,
  !   including a header and a summary line, are written
  ! * At the moment, the following items are not written:
  !   velocity, offshellParameter, coulombPotential, perturbative
  !****************************************************************************
  interface WriteParticle
     module procedure WriteParticle1,WriteParticle2
  end interface


  public :: Write_ReadingInput, Write_InitStatus, writeFileDocu
  public :: WriteParticle, WriteParticle_debug, WriteParticleVector
  public :: intToChar, intToChar4, intToChar_pm
  public :: realToChar, realTochar4
  public :: timeMeasurement, printTime
  public :: DoPr
  public :: notInRelease
  public :: setPrintParticleVectorsFormat

  ! local copy of variable in inputGeneral
  integer, save:: printParticleVectorsFormat = 1


contains

  !****************************************************************************
  !****s* output/notInRelease
  ! NAME
  ! subroutine notInRelease(feature)
  !
  ! PURPOSE
  ! Prints a message which states the given "feature" is not included in the present
  ! released code version and stops the code.
  !
  ! INPUTS
  ! * character(*) :: feature
  !****************************************************************************
  subroutine notInRelease(feature)
    use CallStack, only: TraceBack

    character(*), intent(in) :: feature
    write(*,*)
    write(*,*) line
    write(*,'(3A)') 'Error: The feature ',trim(feature),' is not yet included in this release.'
    write(*,'(A)')  'This feature needs further testing and first results are not yet published.'
    write(*,'(A)')  'Please stay tuned for the next release. Thanks, the GiBUU team !'
    write(*,*)
    call traceback(user_exit_code=-1)
    write(*,*)
    write(*,*) line
    write(*,*)
    stop
  end subroutine notInRelease


  !****************************************************************************
  !****s* output/DoPr
  ! NAME
  ! function DoPr(iLevel)
  ! PURPOSE
  ! Return the Flag, whether a message of the given Error-Level (cf. DoPrLevel)
  ! should be written or not
  ! INPUTS
  ! * integer :: iLevel -- Print-Level
  ! OUTPUT
  ! * logical :: DoPr
  ! NOTES
  ! checks, whether the given iLevel is a valid array-index are omitted
  ! (compile your code with the approbiate flags to ensure array-bound checks)
  !****************************************************************************
  logical function DoPr(iLevel)

    integer, intent(in) :: iLevel

!    DoPr = .FALSE.            ! may be ommitted if you compile with "checks"
!    if (iLevel.lt.-10) return ! may be ommitted if you compile with "checks"
!    if (iLevel.gt.5)   return ! may be ommitted if you compile with "checks"

    DoPr = DoPrLevel(iLevel)

  end function DoPr


  !****************************************************************************
  !****s* output/intToChar
  ! NAME
  ! function intToChar(nr)
  ! PURPOSE
  ! Convert an integer into a character string of length 3.
  ! INPUTS
  ! * integer :: nr -- number to convert (0..999)
  ! OUTPUT
  ! * character(3)
  !****************************************************************************
  character(3) function intToChar (nr)
    integer, intent(in) :: nr

    integer p0,p1,p2
    p2 = Mod(nr/100,10)
    p1 = Mod(nr/10,10)
    p0 = Mod(nr,10)

    intToChar = Achar(p2+48)//Achar(p1+48)//Achar(p0+48)
  end function intToChar

  !****************************************************************************
  !****s* output/intToChar4
  ! NAME
  ! function intToChar4(nr)
  ! PURPOSE
  ! Convert an integer into a character string of length 4.
  ! INPUTS
  ! * integer :: nr -- number to convert (0..9999)
  ! OUTPUT
  ! * character(3)
  !****************************************************************************
  character(4) function intToChar4 (nr)
    integer, intent(in) :: nr

    integer p0,p1,p2,p3
    p3 = Mod(nr/1000,10)
    p2 = Mod(nr/100,10)
    p1 = Mod(nr/10,10)
    p0 = Mod(nr,10)

    intToChar4 = Achar(p3+48)//Achar(p2+48)//Achar(p1+48)//Achar(p0+48)
  end function intToChar4


  !****************************************************************************
  !****s* output/intToChar_pm
  ! NAME
  ! function intToChar_pm(nr)
  ! PURPOSE
  ! Convert an integer into a signed character string of length 2.
  ! INPUTS
  ! * integer :: nr -- number to convert (-9..9)
  ! OUTPUT
  ! * character(3)
  !****************************************************************************
  character(2) function intToChar_pm (nr)
    integer, intent(in) :: nr

    integer p0
    p0=abs(Mod(nr,10))

    if (nr<0) then
       intToChar_pm='-'//Achar(p0+48)
    else
       intToChar_pm='+'//Achar(p0+48)
    end if
  end function intToChar_pm


  !****************************************************************************
  !****s* output/realToChar
  ! NAME
  ! function realToChar(nr)
  ! PURPOSE
  ! Convert a real number into a character string of length 3.
  ! INPUTS
  ! * real :: nr -- number to convert (0..999)
  ! OUTPUT
  ! * character(3)
  !****************************************************************************
  character(3) function realToChar (nr)
    real, intent(in) :: nr

    integer p0,p1,p2
    p2 = int(nr/100)
    p1 = int(nr/10)-10*p2
    p0 = int(nr)-10*p1-100*p2

    realToChar = Achar(p2+48)//Achar(p1+48)//Achar(p0+48)
  end function realToChar


  !****************************************************************************
  !****s* output/realToChar4
  ! NAME
  ! function realToChar4(nr)
  ! PURPOSE
  ! Convert a real number into a character string of length 4.
  ! INPUTS
  ! * real :: nr -- number to convert (0..9999)
  ! OUTPUT
  ! * character(4)
  !****************************************************************************
  character(4) function realToChar4 (nr)
    real, intent(in) :: nr

    integer p3
    p3 = int(nr/1000)

    realToChar4 = Achar(p3+48)//realToChar(nr-1000.*real(p3))
  end function realToChar4


  !****************************************************************************
  !****s* output/Write_ReadingInput
  ! NAME
  ! subroutine Write_ReadingInput(Text,Code,ios)
  ! PURPOSE
  ! Write a Message about starting/finishing "reading" something.
  ! If parameter "ios" is given, this is interpreted as the IOSTAT value
  ! of the reading of the namelist (here parameter "Code" is ignored).
  ! INPUTS
  ! * character(*)     :: Text -- Name of file/Jobcard/...
  ! * integer          :: Code -- 0: start reading, !=0 : reading finished
  ! * integer,optional :: ios  -- error code of reading (i.e. IOSTAT value)
  !****************************************************************************
  subroutine Write_ReadingInput(Text,Code,ios)

    character(*), intent(in) :: Text
    integer, intent(in)      :: Code
    integer, intent(in), optional :: ios

    if (present(ios)) then
       if (ios > 0) then
          write(*,1002) Text
          stop
       else if (ios < 0) then
          write(*,1003) Text
       end if
    else
       if (Code.eq.0) then
          write(*,1000) Text
       else
          write(*,1001) Text
       end if
    end if

1000 FORMAT (/,'------ Init "',A,'": reading...')
1001 FORMAT ('------ Init "',A,'": reading finished.',/)
1002 FORMAT (/,'--- !!!!! ERROR while reading namelist "',A,'" !!!!! STOPPING !!',/)
1003 FORMAT ('--- *** namelist "',A,'" missing in jobcard. Values-->Default!')

  end subroutine Write_ReadingInput


  !****************************************************************************
  !****s* output/Write_InitStatus
  ! NAME
  ! subroutine Write_InitStatus(Text,Code)
  ! PURPOSE
  ! Write a Message about starting/finishing "initializing" something
  ! INPUTS
  ! * character(*) :: Text -- Name of file/Jobcard/...
  ! * integer      :: Code -- 0   : start initialization,
  !                           <>0 : finished initialization
  !****************************************************************************
  subroutine Write_InitStatus(Text,Code)

    character(*), intent(in) :: Text
    integer, intent(in)      :: Code

    if (Code.eq.0) then
       write(*,1000) Text
    else
       write(*,1001) Text
    end if

1000 FORMAT (/,'------ Init "',A,'": ...')
1001 FORMAT ('------ Init "',A,'": finished.',/)

  end subroutine Write_InitStatus



  !****************************************************************************
  ! cf. interface WriteParticle1
  !****************************************************************************
  subroutine WriteParticle1(iFile,iEnsemble,iPart, Part)
    use particleDefinition
    use history, only: history_getParents

    integer,intent(in)          :: iFile
    integer,intent(in),optional :: iEnsemble,iPart
    type(particle),intent(in),optional :: Part
    !character(140), parameter :: FormatPart1 = &
    !   & '(i3,"|",i5,"[",i9,"]:",i3,i3,l3,"  ",1P,e11.3," ",4e11.3,"  ",3e11.3,"  ",0P,f5.3," ",i9,"[",2i9,"] ",1P,e11.3)'

!    character*(*), parameter :: FormatMom = '4e11.3'
    character*(*), parameter :: FormatMom = '4e15.7'

    character(250), parameter :: FormatPart2 = &
       & '(i6,"|",i6,"[",i9,"]:",i3,i3,l2," ",f7.3,1P," |",'//FormatMom// &
       & ',"|",3e11.3,"| ",0P,f5.3," ",i9,"[",2i9,"] ",'&
       & //'1P,e11.3," time:",0P,3f9.3,l2," parents:",3i4)'

    integer :: parents(1:3)

!!$    write(iFile,FormatPart1) iEnsemble,iPart, &
!!$         Part%number, Part%ID, Part%charge, Part%antiparticle, &
!!$         Part%mass, Part%momentum, Part%position,  &
!!$         Part%ScaleCS, Part%firstEvent, Part%event, &
!!$         Part%perWeight

    if (present(iEnsemble) .and. present(iPart) .and. present(Part)) then
       parents = history_getParents(Part%history)
       write(iFile,FormatPart2) iEnsemble,iPart, &
            Part%number, Part%ID, Part%charge, Part%antiparticle, &
            Part%mass, Part%momentum, Part%position,  &
            min(9.999,Part%ScaleCS), Part%firstEvent, Part%event, &
            Part%perWeight, &
            max(-9.9E0,Part%formationTime), &
            Part%productionTime, Part%lastCollisionTime, &
            Part%in_Formation, parents
    else
       ! Just write the header
       write(iFile,'(2A)') 'Ensemb| Index[   Number]: ID IQ Anti mass | momentum(0:3) | position(1:3) | ', &
       & 'ScaleCS firstEvent [event] perWeight, max(-9.9d0,formationTime), productionTime,lastCollisionTime, in_Formation, parents'
    end if

  end subroutine WriteParticle1
  !-------------------------------------------------------------------------
  subroutine WriteParticle2(iFile,iEnsemble,Part)
    use particleDefinition

    integer,intent(in)        :: iFile
    integer,intent(in)        :: iEnsemble
    type(particle),dimension(:),intent(in) :: Part

    type(particle) :: hPart

    integer :: i

    call setToDefault(hPart)

    write(iFile,'(200("="))')
    call WriteParticle1(iFile)

    do i=1,size(Part)
       if (Part(i)%ID<=0) cycle

       hPart%charge  = hPart%charge  + Part(i)%charge
       hPart%momentum= hPart%momentum+ Part(i)%momentum

       call WriteParticle1(iFile,iEnsemble,i,Part(i))
    end do

    write(iFile,'(200("-"))')

    hPart%mass = sqrtS(hPart)
    call WriteParticle1(iFile,iEnsemble,0,hPart)
    write(iFile,'(200("="))')

  end subroutine WriteParticle2


  !****************************************************************************
  !****s* output/WriteParticle_debug
  ! NAME
  ! subroutine WriteParticle_debug (Part, med)
  ! PURPOSE
  ! Write the stored information about a particle to a file/stdout.
  ! INPUTS
  ! * type(particle) :: Part          --- Particle to print
  ! * type(medium), optional :: med   --- medium information (optional)
  !****************************************************************************
  subroutine WriteParticle_debug (Part, med)
    use particleDefinition
    use mediumDefinition
    use minkowski, only: abs4

    type(particle),intent(in) :: Part
    type(medium), intent(in), optional :: med

    write(*,'(79("-"))')
    write(*,*) '#:', Part%number, 'ID:',Part%ID, 'Charge:',Part%charge, 'Antiparticle:',Part%antiparticle
    write(*,*) 'Perturbative?', part%perturbative
    write(*,*) 'bare Mass=', Part%mass
    write(*,*) 'eff. Mass=', abs4(Part%momentum)

    write(*,'(A,1P,4E17.9)') ' Momentum=', Part%momentum
    write(*,'(A,1P,17(" "),3E17.9)') ' Position=', Part%position
    write(*,'(A,1P,17(" "),3E17.9)') ' Velocity=', Part%velocity

    write(*,*) 'Offshellparameter=', Part%offshellparameter
    write(*,*) 'History=',Part%history
    write(*,*) 'ScaleCS=',Part%ScaleCS, 'FirstEvent=',Part%firstEvent
    write(*,*) 'Event(1:2)=', Part%event,  'Perweight=',Part%perWeight

    if (present(med)) write(*,*) 'Density=', med%density
    write(*,'(79("-"))')
    write(*,*)

  end subroutine WriteParticle_debug


  !****************************************************************************
  !****s* output/writeFileDocu
  ! NAME
  ! subroutine writeFileDocu(filename,docu)
  !
  ! PURPOSE
  ! Writes the information of a given filename and its documenation to a central file.
  !
  ! INPUTS
  ! * character(*) :: filename
  ! * character(*) :: docu -- Documentation of the output file "filename"
  !****************************************************************************
  subroutine writeFileDocu(filename,docu)

    character(*) :: filename
    character(*) :: docu       ! documentation
    logical, save :: firstTime=.true.

    if (firstTime) then
       open(20,file="README.OutputFiles.txt")
       firstTime=.false.
    else
       open(20,file="README.OutputFiles.txt",position='Append')
    end if
    write(20,*)
    write(20,'(2A)') 'File name: ', filename
    write(20,'(2A)') 'Purpose:   ', docu
    write(20,*)
    write(20,*) line
    close(20)

  end subroutine writeFileDocu


  !****************************************************************************
  !****s* output/WriteParticleVector
  ! NAME
  ! subroutine WriteParticleVector(filename,pv)
  !
  ! PURPOSE
  ! Write a particle-vector (e.g. PertParticles) to files
  !
  ! Opens 7 files and writes the information of the particles to these files.
  ! Only for the file "...ALL.dat" the subroutine WriteParticle is used; all
  ! other information is written in an own format.
  !
  ! INPUTS
  ! * character(*) :: filename -- base name of files (some extensions are added)
  ! * type(particle),dimension(:,:) :: pv -- particle vector
  !****************************************************************************
  subroutine WriteParticleVector(filename,pv)
    use particleDefinition
    use CallStack, only: TraceBack

    character(*), intent(in) :: filename
    type(particle),dimension(:,:),intent(in) :: PV ! particleVector

    select case (printParticleVectorsFormat)
    case (1)
       call doASCII
    case (2)
       call doBinary
    case default
       call TraceBack("wrong printParticleVectorsFormat")
    end select

  contains
    subroutine doASCII

      integer :: i,j

      open(96,File=filename//"_velo.dat")
      open(97,File=filename//"_pos.dat")
      open(98,File=filename//"_mom.dat")
      open(99,File=filename//"_1.dat")
      open(100,File=filename//"_forma.dat")
      open(101,File=filename//"_2.dat")
      open(102,File=filename//"_ALL.dat")

      write(96,'(A)') "# number,ensemble, ID of particle,charge, velocity "
      write(97,'(A)') "# number,ensemble, ID of particle,charge, position[fm]"
      write(98,'(A)') "# number,ensemble, ID of particle,charge, 4-momentum[GeV]"
      write(99,'(A)') "# number,ensemble, event, perturbative, antiparticle,perweight"
      write(100,'(A)')"# number,ensemble, scaleCS, ,formationTime, productionTime,lastCollisionTime, in_Formation"
      write(101,'(A)')"# number,ensemble, mass,coulombPotential, offShellParameter, firstevent "
      call WriteParticle(102)

      do i=1,size(pv,dim=1)
         do j=1,size(pv,dim=2)
            if (pv(i,j)%Id <= 0) cycle
            write(96,'(4I10,3E14.5)')    pv(i,j)%number, i, pv(i,j)%ID,pv(i,j)%Charge,pv(i,j)%velocity
            write(97,'(4I10,3E14.5)')    pv(i,j)%number, i, pv(i,j)%ID,pv(i,j)%Charge,pv(i,j)%position
            write(98,'(4I10,4E14.5)')    pv(i,j)%number, i, pv(i,j)%ID,pv(i,j)%Charge,pv(i,j)%momentum

            write(99,'(4I10,2L4,e14.3)')  pv(i,j)%number, i, pv(i,j)%event,pv(i,j)%perturbative,&
                 & pv(i,j)%antiparticle,pv(i,j)%perweight
            write(100,'(2I10,4F10.4,L4)') pv(i,j)%number, i, pv(i,j)%scaleCS,pv(i,j)%formationTime,pv(i,j)%productionTime,&
                 & pv(i,j)%lastCollisionTime,pv(i,j)%in_Formation
            write(101,'(2I10,3F10.4,I8)') pv(i,j)%number, i, pv(i,j)%mass,0.0,pv(i,j)%offShellParameter,&
                 & pv(i,j)%firstEvent

            call WriteParticle(102,i,j,pv(i,j))
         end do
      end do

      close(96)
      close(97)
      close(98)
      close(99)
      close(100)
      close(101)
      close(102)

    end subroutine doASCII

    subroutine doBinary

      integer:: i,j,n

      open(96,file=filename//".bin",status='UNKNOWN',form='UNFORMATTED')
      rewind(96)
      write(96) size(pv,dim=1)

      do i=1,size(pv,dim=1)
         n = size(pv,dim=2)
         do j=n,1,-1
            if (pv(i,j)%Id <= 0) cycle
            n = j
            exit
         end do

         write(96) n
         do j=1,n
            write(96) pv(i,j)
         end do

      end do
      close(96)

    end subroutine doBinary

  end subroutine WriteParticleVector

  !****************************************************************************
  !****s* output/setPrintParticleVectorsFormat
  ! NAME
  ! subroutine setPrintParticleVectorsFormat(f)
  !
  ! PURPOSE
  ! set the format for printParticleVector
  !****************************************************************************
  subroutine setPrintParticleVectorsFormat(f)
    integer, intent(in) :: f
    printParticleVectorsFormat = f
  end subroutine setPrintParticleVectorsFormat


  !****************************************************************************
  !****s* output/timeMeasurement
  ! NAME
  ! subroutine timeMeasurement(ForceReset)
  !
  ! PURPOSE
  ! Use this as an stopwatch.
  !
  ! If you do not use the optional parameter "ForceReset", the default
  ! behaviour is like this:
  ! * the first call sets the stopwatch to zero
  ! * all following (second, third, ...) calls just print the time
  !
  ! only by calling this routine and setting the optional parameter
  ! ForceReset=.true. can rewind the stopwatch to zero
  !
  ! INPUTS
  ! * logical,optional :: ForceReset -- see above
  ! * integer,optional :: iFile -- File to print the message
  !
  ! OUTPUT
  ! on stdout
  !****************************************************************************
  subroutine timeMeasurement(ForceReset,iFile)

    logical, intent(in), optional :: ForceReset
    integer, intent(in), optional :: iFile

    integer :: timeEnd, timeRate
    logical, save :: DoReset=.true. ! Initialize at first time
    integer, save :: timestart
    integer :: ii

    if (present(ForceReset)) then
      DoReset=ForceReset
    end if
    ii = 6
    if (present(iFile)) ii = iFile

    if (DoReset) then
      call system_clock(timeEnd,timeRate)
      timeStart = timeEnd
      DoReset = .false.
    else
      call system_clock(timeEnd,timeRate)   ! time measurement for performance
      write(ii,'(1X,"#",10(" "),"Time elapsed : ",f15.4," seconds")') real(timeEnd-timeStart)/real(timeRate)
    end if

  end subroutine timeMeasurement


  !****************************************************************************
  !****s* output/PrintTime
  ! NAME
  ! subroutine PrintTime
  ! PURPOSE
  ! print out date and time information
  !****************************************************************************
  subroutine PrintTime(c)
    character(len=*),intent(in)::c
    integer:: v(8)
    call date_and_time(values=v)
    write(*,'(A,": ",i2.2,".",i2.2,".",i4.4,"  ",i2.2,":",i2.2,":",i2.2)') c,v(3),v(2),v(1),v(5),v(6),v(7)
  end subroutine PrintTime


end module output
