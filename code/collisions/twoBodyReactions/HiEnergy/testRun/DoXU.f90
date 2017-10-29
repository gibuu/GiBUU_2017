program DoXU

  use output
  use version
  use particleDefinition
  use particleProperties
  use hadronFormation
  use CollTools
  use VMMassPythia

  IMPLICIT NONE

  integer :: iEV
!  integer,parameter :: NEV=10000000
  integer,parameter :: NEV=10

  COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
  integer N,NPAD,K
  double precision P,V
  SAVE /PYJETS/

  integer nArrMax
  parameter (nArrMax=200)   ! maximum size of arrays

  common /DataGJV/&
       &     Arr(3,4,nArrMax),&    ! 3* 4D-Vertizes
       &     EArr(6,nArrMax),&     ! errFlag, rank
       &     verb,&                ! verbosity
       &     AtOrigin             ! treatment of outmost prod points

  double precision Arr
  integer EArr
  integer verb
  logical AtOrigin
  
  save /DataGJV/

  call PrintVersion

  call forceInitFormation
  call InitParticleProperties


  !--------------------
  !... Initializations:
  !--------------------

  call SetSomeDefaults_PY

  call Init_VM_Mass(200d0)

  if (.not.useJetSetVec) then
     write(*,*) 'set useJetSetVec! stop.'
     stop
  end if

  CALL PYINIT('CMS','p','n', 200d0)

  verb=1

  do iEV=1,NEV
     call GetJetsetVecINIT
     CALL PYEVNT
     verb=1
     call GetJetsetVec(.TRUE.)

     call PYLIST(2)
     call GetJetSetVec_List(6,1,N)

     write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

     call SFREP_Write(1,N)
     write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'


     call SFREPS_Write(0)
     write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     call SFREPS_Write(1)
     write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     call SFREPS_Write(2)
     write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     call SFREPS_Write(3)
     

     stop

     call GetJetsetVecCheckT(-1d-5)

     call GetJetsetVecPYEDIT

     
     call PYLIST(2)
     stop
  end do

end program DoXU
