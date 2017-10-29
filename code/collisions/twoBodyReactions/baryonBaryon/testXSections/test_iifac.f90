! This is a test case to compare the tabulatied BB->BB isospin factors
! (iifac0/iifac1) to a direct calculation (isf), and to do some consistency checks.

program test_iifac

  use IDtable, only: nucleon, Delta
  use inputGeneral, only: readInputGeneral
  use particleProperties, only: hadron, initParticleProperties
  use barBar_BarBar, only: iifac0, iifac1

  implicit none

  integer :: IDout(1:2),is(1:2), q1, q2
  real :: sum1, sum2

  call readInputGeneral()
  call initParticleProperties()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,"**********************************************************************"
  print *,"N N -> N N"
  is = (/1,1/)
  IDout = (/Nucleon,Nucleon/)

  call compare(is,IDout,(/1,1/),"p p -> N N")
  call compare(is,IDout,(/1,0/),"p n -> N N")
  call compare(is,IDout,(/0,0/),"n n -> N N")

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,"**********************************************************************"
  print *,"N N -> N D"
  is = (/1,1/)
  IDout = (/Nucleon,Delta/)

  call compare(is,IDout,(/1,1/),"p p -> N D")
  call compare(is,IDout,(/1,0/),"p n -> N D")
  call compare(is,IDout,(/0,0/),"n n -> N D")

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,"**********************************************************************"
  print *,"N N -> D D"
  is = (/1,1/)
  IDout = (/Delta,Delta/)

  call compare(is,IDout,(/1,1/),"p p -> D D")
  call compare(is,IDout,(/1,0/),"p n -> D D")
  call compare(is,IDout,(/0,0/),"n n -> D D")

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,"**********************************************************************"
  print *,"N D -> N N"
  is = (/1,3/)
  IDout = (/Nucleon,Nucleon/)

  call compare(is,IDout,(/1,2/),"p D++ -> N N")
  call compare(is,IDout,(/1,1/),"p D+ -> N N")
  call compare(is,IDout,(/1,0/),"p D0 -> N N")
  call compare(is,IDout,(/1,-1/),"p D- -> N N")
  call compare(is,IDout,(/0,2/),"n D++ -> N N")
  call compare(is,IDout,(/0,1/),"n D+ -> N N")
  call compare(is,IDout,(/0,0/),"n D0 -> N N")
  call compare(is,IDout,(/0,-1/),"n D- -> N N")

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,"**********************************************************************"
  print *,"N D -> N D"
  is = (/1,3/)
  IDout = (/Nucleon,Delta/)

  call compare(is,IDout,(/1,2/),"p D++ -> N D")
  call compare(is,IDout,(/1,1/),"p D+ -> N D")
  call compare(is,IDout,(/1,0/),"p D0 -> N D")
  call compare(is,IDout,(/1,-1/),"p D- -> N D")
  call compare(is,IDout,(/0,2/),"n D++ -> N D")
  call compare(is,IDout,(/0,1/),"n D+ -> N D")
  call compare(is,IDout,(/0,0/),"n D0 -> N D")
  call compare(is,IDout,(/0,-1/),"n D- -> N D")

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,"**********************************************************************"
  print *,"D D -> N N"
  is = (/3,3/)
  IDout = (/Nucleon,Nucleon/)

  call compare(is,IDout,(/2,2/),"D++ D++ -> N N")
  call compare(is,IDout,(/2,1/),"D++ D+ -> N N")
  call compare(is,IDout,(/2,0/),"D++ D0 -> N N")
  call compare(is,IDout,(/2,-1/),"D++ D- -> N N")
  call compare(is,IDout,(/1,1/),"D+ D+ -> N N")
  call compare(is,IDout,(/1,0/),"D+ D0 -> N N")
  call compare(is,IDout,(/1,-1/),"D+ D- -> N N")
  call compare(is,IDout,(/0,0/),"D0 D0 -> N N")
  call compare(is,IDout,(/0,-1/),"D0 D- -> N N")
  call compare(is,IDout,(/-1,-1/),"D- D- -> N N")

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,"**********************************************************************"

contains


  subroutine compare(is,IDout,q,str)
    integer, intent(in), dimension(1:2) :: is, q, IDout
    character(*), intent(in) :: str
    integer :: q1min,q1max,q2min,q2max

    ! Charges which are possible for the first final state particle
    q1min=int(float(1-hadron(idOut(1))%isoSpinTimes2)/2.)
    q1max=int(float(1+hadron(idOut(1))%isoSpinTimes2)/2.)
    ! Charges which are possible for the second final state particle
    q2min=int(float(1-hadron(idOut(2))%isoSpinTimes2)/2.)
    q2max=int(float(1+hadron(idOut(2))%isoSpinTimes2)/2.)

    print *,"************************************"
    print *,str
    sum1 = iifac0(IDout,is,q)
    sum2 = isf(is,q,IDout)
    print '(2F7.4,L2)', sum1, sum2, abs(sum1-sum2)<1E-6
    print *,"------"
    do q1=q1min,q1max
      do q2=q2min,q2max
         call compare1(is,IDout,q,(/q1,q2/))
      end do
    end do
    print *,"------"
    print '(2L2)', abs(sum1)<1E-6, abs(sum2)<1E-6

  end subroutine


  subroutine compare1(is, IDout, q, q_out)
    integer, intent(in), dimension(1:2) :: is, q, IDout,q_out
    real :: i1, i2
    if (sum(q_out)/=sum(q)) return
    i1 = iifac1(IDout,is,q,q_out(1),q_out(2))
    i2 = isf(is,q,IDout,q_out)
    print '(2F7.4,L2)', i1, i2, abs(i1-i2)<1E-6
    sum1=sum1-i1
    sum2=sum2-i2
  end subroutine


  !*****************************************************************************
  !****f* test_iifac/isf
  ! NAME
  ! real function isf (is_in, q_in, IDout, q_out)
  ! NOTES
  ! This function calculates the isospin factors given by equ. (A.5) of
  ! Effenberger's PhD thesis, times a factor for identical final state particles.
  ! If "q_out" is not given, we sum over all possible final state charges.
  ! This routine is only valid for B_1 B_2 -> B_3 B_4,
  ! where all B's are non-strange/non-charmed Baryons.
  ! INPUTS
  ! * is_in  : isospins of initial state times 2
  ! * q_in   : charges of initial state particles
  ! * IDout  : IDs of final state
  ! * q_out  : charges of final state particles (optional)
  !*****************************************************************************
  real function isf (is_in, q_in, IDout, q_out)

    use particleProperties, only: hadron, isBaryon
    use clebschGordan, only: CG

    integer, intent(in), dimension(1:2)           :: is_in, q_in, IDout
    integer, intent(in), dimension(1:2), optional :: q_out
    integer, dimension(2) :: is_out, iz_in, iz_out
    integer :: M, I, iz_out1, iz_out2

    isf = 0.

    is_out(1) = hadron(IDout(1))%isoSpinTimes2
    is_out(2) = hadron(IDout(2))%isoSpinTimes2

    iz_in  = 2*q_in - 1   ! assuming non-strange baryons
    M = sum(iz_in)

    if (present(q_out)) then

      ! final state charges given as input
      if (sum(q_in) /= sum(q_out)) return
      iz_out = 2*q_out - 1   ! assuming non-strange baryons

      do I = abs(M),min(sum(is_in),sum(is_out))
         isf = isf + CG ( is_in(1), is_in(2),I, iz_in(1), iz_in(2),M)**2 * &
              CG (is_out(1),is_out(2),I,iz_out(1),iz_out(2),M)**2
      end do

    else

      ! sum over final state charges
      do I = abs(M),min(sum(is_in),sum(is_out))
         do iz_out1 = -is_out(1),is_out(1),2  ! 2-steps
            iz_out2 = M-iz_out1
            isf = isf + CG ( is_in(1), is_in(2),I,iz_in(1),iz_in(2),M)**2 * &
                 CG (is_out(1),is_out(2),I, iz_out1, iz_out2,M)**2
         end do
      end do

    end if

    ! factor for identical particle IDs in final state
    if (IDout(1)==IDout(2)) isf = isf * 0.5

  end function isf


end program
