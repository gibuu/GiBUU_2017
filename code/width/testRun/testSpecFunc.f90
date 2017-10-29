program test
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: InitParticleProperties
  implicit none

  call readinputGeneral
  call InitParticleProperties
  call tester

contains

  subroutine tester
    use spectralFunc, only: specFunc
    use IDTABLE

    integer, parameter :: ID=Delta
    integer, parameter :: charge=1
    real, dimension(0:3) :: p,r
    real :: m, pAbs,intgrl,dm,sf
    real, parameter :: m_min=0.
    real, parameter :: m_max=6
    integer :: numSteps_m,i,j

    r = 0.
    dm=0.0005
    numSteps_m=NINT((m_max-m_min)/dm)

    do i=0,50
      pAbs=float(i)*0.02
      p(1:3) = (/ 0., 0., pAbs /)
      intgrl=0
      write(*,*) i
      do j=0,numSteps_m
        m=0.+float(j)*dm
        p(0)=sqrt(m**2+pabs**2)
        sf = specFunc(ID,charge,p,r)
        write(14,*) m, pabs, sf
        intgrl = intgrl + sf*dm*2*m
      end do
      write(14,*)
      write(12,*) pAbs,intgrl
    end do
  end subroutine tester

end program
