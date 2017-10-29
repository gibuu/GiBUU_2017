!******************************************************************************
!****m* /NuclearPDF
! NAME
! module NuclearPDF
! PURPOSE
! Implement nuclear parton distributions
!******************************************************************************
module NuclearPDF

  implicit none
  private

  public :: UseNuclearPDF
  public :: SetNuclearPDFA
  public :: DoNuclearPDF

  !****************************************************************************
  !****s* NuclearPDF/NuclearPDFtype
  ! SOURCE
  !
  integer,  save :: NuclearPDFtype = 0
  !
  ! PURPOSE
  ! Select, which nuclear modification of the parton distribution functions
  ! is used:
  ! * 0: no modification
  ! * 1: EKS 98
  !****************************************************************************

  logical, save :: initFlag = .true.

  integer, save :: NuclearPDFA = 1

contains

  !****************************************************************************
  !****s* NuclearPDF/readInput
  ! NAME
  ! subroutine readInput
  !
  ! PURPOSE
  ! Reads input in jobcard out of namelist "NuclearPDF"
  !****************************************************************************
  subroutine readInput
    use output
    use CallStack

    integer :: ios

    NAMELIST /NuclearPDF/ NuclearPDFtype

    if (.not.initFlag) return

    call Write_ReadingInput('NuclearPDF',0)
    rewind(5)
    read(5,nml=NuclearPDF,IOSTAT=ios)
    call Write_ReadingInput("NuclearPDF",0,ios)

    select case (NuclearPDFtype)
    case (0)
       write(*,*) 'Nuclear PDF: none'
    case (1)
       write(*,*) 'Nuclear PDF: EKS98'
    case default
       call TRACEBACK('wrong NuclearPDFtype')
    end select

    call Write_ReadingInput('NuclearPDF',1)
    initFlag = .false.

  end subroutine readInput

  !****************************************************************************
  !****f* NuclearPDF/UseNuclearPDF
  ! NAME
  ! logical function UseNuclearPDF()
  !
  ! PURPOSE
  ! Return true, if we will modify the pdf for the nucleus
  !****************************************************************************
  logical function UseNuclearPDF()
    if (initFlag) call readInput
    UseNuclearPDF = (NuclearPDFtype .ne. 0).and.(NuclearPDFA .gt. 1)
  end function UseNuclearPDF

  !****************************************************************************
  !****s* NuclearPDF/SetNuclearPDFA
  ! NAME
  ! subroutine SetNuclearPDFA(A)
  !
  ! PURPOSE
  ! Set the nucleus to be used
  !
  ! INPUTS
  ! * integer :: A -- nuclear size to be used
  ! OUTPUT
  ! none
  !****************************************************************************
  subroutine SetNuclearPDFA(A)
    integer, intent(in) :: A
    NuclearPDFA = A
  end subroutine SetNuclearPDFA

  !****************************************************************************
  !****s* NuclearPDF/DoNuclearPDF
  ! NAME
  ! subroutine DoNuclearPDF(x,Q, xPQ)
  !
  ! PURPOSE
  ! Apply the nuclear modifications.
  !
  ! INPUTS
  ! * real:: x,Q -- arguments of pdf
  ! * real, dimension(-25:25) :: proton pdfs as in PYPDFU
  ! OUTPUT
  ! * real, dimension(-25:25) :: modified pdfs
  !****************************************************************************
  subroutine DoNuclearPDF(x,Q, xPQ)
    double precision, intent(in) :: X, Q
    double precision, dimension(-25:25), intent(inOut) :: xPQ

    real :: Ruv,Rdv,Rub,Rdb,Rs,Rc,Rb,Rt,Rg

    if (.not.UseNuclearPDF()) return

    select case (NuclearPDFtype)
    case (1)
       call eks98(x,Q,real(NuclearPDFA),Ruv,Rdv,Rub,Rdb,Rs,Rc,Rb,Rt,Rg)

       XPQ(0) = XPQ(0)*Rg
       XPQ(1) = (XPQ(1)-XPQ(-1))*Rdv+XPQ(-1)*Rdb ! = DNV+DSEA
       XPQ(-1)= XPQ(-1)*Rdb                      ! = DSEA
       XPQ(2) = (XPQ(2)-XPQ(-2))*Ruv+XPQ(-2)*Rub ! = UPV+USEA
       XPQ(-2)= XPQ(-2)*Rub                      ! = USEA
       XPQ(3) = XPQ(3)*Rs
       XPQ(-3)= XPQ(-3)*Rs
       XPQ(4) = XPQ(4)*Rc
       XPQ(-4)= XPQ(-4)*Rc
       XPQ(5) = XPQ(5)*Rb
       XPQ(-5)= XPQ(-5)*Rb
       XPQ(6) = XPQ(6)*Rt
       XPQ(-6)= XPQ(-6)*Rt

    end select


  end subroutine DoNuclearPDF



end module NuclearPDF
