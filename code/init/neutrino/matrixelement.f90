!******************************************************************************
!****m* /NeutrinoMatrixElement
! NAME
! module NeutrinoMatrixElement
!
! PURPOSE
! This module calculates the squared spin summed and averaged
! matrix element for neutrino induced reactions.
!
!******************************************************************************

module NeutrinoMatrixElement

  implicit none
  private

  !****************************************************************************
  !****g* NeutrinoMatrixElement/which_resonanceModel
  ! SOURCE
  integer, save :: which_resonanceModel=0
  ! PURPOSE
  ! to change between different realizations of the matrix elements:
  ! * 0 = with Fortran calculated matrix elements containing all resonances (default)
  ! * 1 = with Mathematica calculated matrix elements (only Delta)
  ! * 2 = Rein and Sehgals matrix elements
  !****************************************************************************

  logical, save :: initFlag=.true.


  public :: nuMaEl



contains

  subroutine readInput_nu_MatrixElement
    use output

    integer :: ios

    !**************************************************************************
    !****n* NeutrinoMatrixElement/neutrino_matrixelement
    ! NAME
    ! NAMELIST neutrino_matrixelement
    ! PURPOSE
    ! Includes parameters for neutrino matrix elements:
    ! * which_resonanceModel
    !**************************************************************************

    NAMELIST /neutrino_matrixelement/  which_resonanceModel
    call Write_ReadingInput('neutrino_matrixelement',0)
    rewind(5)
    read(5,nml=neutrino_matrixelement,IOSTAT=ios)
    call Write_ReadingInput("neutrino_matrixelement",0,ios)
    call Write_ReadingInput('neutrino_matrixelement',1)

    if (which_resonanceModel.eq.0) then
       write(*,*) 'default choice of matrixelements'
    else if (which_resonanceModel.eq.1) then
       write(*,*) 'old matrixelement calculated by Mathematica'
    else if (which_resonanceModel.eq.2) then
       write(*,*) 'Rein and Sehgal matrixelements'
    end if

  end subroutine readInput_nu_MatrixElement


  !****************************************************************************
  !****f* NeutrinoMatrixElement/nuMaEl
  ! NAME
  ! real function nuMaEl(process_ID,finalstate_ID,charge_in,k_in,k_out,p_in,p_out,bareMass,position)
  ! PURPOSE
  ! calculate the neutrino matrix element
  ! INPUTS
  ! ...
  ! OUTPUT
  ! ...
  !****************************************************************************
  real function nuMaEl(process_ID,finalstate_ID,charge_in,k_in,k_out,p_in,p_out,bareMass,position)
    use idtable
    use ParticleProperties, only: hadron
    use minkowski, only: contract, SP
    !use neutrino_IDTable
    use hadronTensor_ResProd
    use matrixElementQE
    use nu_MatrixElementRS
    use leptonTensor
    use leptonicID

    real, dimension(0:3), intent(in) :: k_in, k_out, p_in, p_out
    integer, intent(in) :: finalstate_ID, charge_in
    integer, intent(in) :: process_ID
    real   , intent(in) :: bareMass            ! bare mass of the produced particle (i.e. without potentials!)
    real, dimension(1:3), intent(in), optional :: position
    complex,dimension(0:3,0:3) :: hadronTensor,leptonTens
    complex :: MaEl
    !real :: s,t,Mi,Mf
    real :: isofactor


    if (initFlag) then
       call readInput_nu_MatrixElement
       initFlag=.false.
    end if

    nuMaEl=0.



    select case (finalstate_ID)

    case (nucleon)
       if (charge_in.eq.proton.and.process_ID.eq.CC) then
          !reaction not possible
          return
       end if
       if (charge_in.eq.neutron.and.process_ID.eq.antiCC) then
          !reaction not possible
          return
       end if
       if (present(position)) then
          nuMaEl=MatrixElementforQE(process_ID,charge_in,k_in,k_out,p_in,p_out,position)
       else
          nuMaEl=MatrixElementforQE(process_ID,charge_in,k_in,k_out,p_in,p_out)
       end if

    case (delta:F37_1950)

       select case (which_resonanceModel)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       case (0) !all resonances based on MAID formfactors
          isofactor=1.
          if (charge_in.eq.proton.and.process_ID.eq.CC) then
             if (hadron(finalstate_ID)%isoSpinTimes2.eq.3) then
                isofactor=3.
             else          !reaction not possible
                return
             end if
          end if
          if (charge_in.eq.neutron.and.process_ID.eq.antiCC) then
             if (hadron(finalstate_ID)%isoSpinTimes2.eq.3) then
                isofactor=3.
             else          !reaction not possible
                return
             end if
          end if

          if (finalstate_ID.eq.S11_2090) return
          if (finalstate_ID.eq.D13_2080) return
          if (finalstate_ID.eq.G17_2190) return
          if (finalstate_ID.eq.P11_2100) return
          if (finalstate_ID.eq.P13_1900) return
          if (finalstate_ID.eq.F15_2000) return
          if (finalstate_ID.eq.S31_1900) return
          if (finalstate_ID.eq.D33_1940) return
          if (finalstate_ID.eq.D35_1930) return
          if (finalstate_ID.eq.D35_2350) return
          if (finalstate_ID.eq.P31_1750) return
          if (finalstate_ID.eq.F35_1750) return

          leptonTens=leptonicTensor(process_ID,k_in,k_out)

          if (hadronTensor_R(p_in,p_out,finalstate_ID,charge_in,process_ID,hadronTensor,bareMass) ) then
             MaEl=Contract( leptonTens, hadronTensor)
          else
             MaEl=0.
          end if

          if (abs(AIMAG(MaEl)).gt.0.0000001) then
             write(*,*) 'ERROR: (Matrix element)^2  is imaginary:', MaEl
          end if
          nuMaEl=isofactor*REAL(MaEl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case (1) !only Delta, uses the matrix elements calculated by Mathematica
          if (finalstate_ID.eq.delta) nuMaEl=nu_MaElDELTA()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case (2) !Rein-Sehgal matrix elements
          nuMaEl=nu_MaEl_RS(process_ID,finalstate_ID,charge_in,k_in,k_out,p_in,p_out)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       case default
          write(*,*) 'which_resonanceModel not implemented',which_resonanceModel
          stop
       end select


    case default
       write(*,*) 'finalstate_ID not yet implemented'
       stop
    end select





  contains




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    real function nu_MaElDELTA()
      use FF_Delta_production
      use constants

      real :: ml
      real :: ampl0,ampl1,ampl2
      real :: c3v,c4v,c5v,c6v,c3a,c4a,c5a,c6a
      real :: c2=0.
      real :: ml2, ml4, Mi2, Mi3, Mi4, Mi5, Mi6, Mf2, Mf3, Mf4, Mf5, Mf6
      real :: Mf7, Mf8, M2, M4, t2, t3, t4, s2
      real :: s,t,Mi,Mf,ml_out,M,Qs
      real :: a
      real :: coupling=0.


      M=mN

      t=SP(k_in-k_out,k_in-k_out)
      Qs=-t
      s=SP(k_in+p_in,k_in+p_in)
      ml_out=sqrt(max(SP(k_out,k_out),0.))
      Mi=sqrt(SP(p_in,p_in))
      Mf=sqrt(SP(p_out,p_out))

      select case (process_ID)
      case (NC)
         coupling=GF**2/4.

      case (CC)
         coupling=GF**2/4.*coscab**2

      case (EM)
         coupling=(4.*pi*alphaQED)**2/Qs**2*1./4./2.
      end select

      select case (process_ID)
      case (NC,CC)
         a=-1.
      case (antiNC,antiCC)
         a=1.
      case (EM,antiEM)
         a=0.
      end select

      call formfactors_Delta(Qs,Mf,process_ID,c3v,c4v,c5v,c6v,c3a,c4a,c5a,c6a)

      ml=ml_out

      ml2=ml**2
      ml4=ml**4
      Mi2=Mi**2
      Mi3=Mi**3
      Mi4=Mi**4
      Mi5=Mi**5
      Mi6=Mi**6
      Mf2=Mf**2
      Mf3=Mf**3
      Mf4=Mf**4
      Mf5=Mf**5
      Mf6=Mf**6
      Mf7=Mf**7
      Mf8=Mf**8
      M2=M**2
      M4=M**4
      t2=t**2
      t3=t**3
      t4=t**4
      s2=s**2


      if (process_ID.eq.CC.and.charge_in.eq.proton) then
         c2=3.
      else if (process_ID.eq.antiCC.and.charge_in.eq.neutron) then
         c2=3.
      else
         c2=1.
      end if

      ampl0=-((1./(3*M2*Mf2))*(8*c2*((ml2 + t)*Mi6 + (ml2 + t)*&
           (2*ml2 + 3*Mf2 - 2*s - 3*t)*Mi4 - &
           (2*(t - 2*Mf2)*ml4 + (Mf4 + 4*s*Mf2 + t2)*ml2 - &
           t*(Mf4 - 8*s*Mf2 - 4*t*Mf2 + 2*s2 + 3*t2 + 4*s*t))*Mi2 - &
           4*(ml4 + t*ml2 - 2*t2)*Mf3*Mi + 2*ml4*Mf2*(t - Mf2) + &
           ml2*(t - Mf2)*(3*Mf4 - 2*t*Mf2 + t2 + 2*s*(t - 3*Mf2)) - &
           t*(-3*Mf6 + 3*t*Mf4 - t2*Mf2 + t3 + 2*s2*(t - 3*Mf2) + &
           2*s*(3*Mf4 - 4*t*Mf2 + t2)))*C3a**2)) + &
           ((16*C6a*ml2*(Mi2 + 2*Mf*Mi + Mf2 - t)*((ml2 + t)*Mi2 + &
           ml2*(t - Mf2) + t*(Mf2 - 2*s - t))*c2)/(3*M**3*Mf) + &
           (1./(3*M*Mf))*(16*C5a*((t - ml2)*Mi4 + 4*t*Mf*Mi3 + &
           2*(-t2 + s*t + ml2*Mf2)*Mi2 + 4*Mf*(ml4 + (Mf2 - s)*ml2 - &
           t*(Mf2 + t))*Mi - 2*ml4*(t - Mf2) + ml2*(-3*Mf2 + 4*s + t)* &
           (t - Mf2) + t*(-3*Mf4 + 2*s*Mf2 + 2*t*Mf2 - 2*s2 + t2 - 2*s*t))* &
           c2) + (1./(3*M**3*Mf))*(8*C4a*((-(ml2 + t))*Mi6 - &
           4*(ml2 + t)*Mf*Mi5 - (ml2 + t)*(2*ml2 + 3*Mf2 - 2*s - 3*t)*Mi4 + &
           8*(-ml4 + s*ml2 + t*(s + t))*Mf*Mi3 + &
           (2*(t - Mf2)*ml4 + (Mf4 + 4*s*Mf2 + 2*t*Mf2 + t2)*ml2 - &
           t*(Mf4 - 8*s*Mf2 + 2*s2 + 3*t2 + 4*s*t))*Mi2 + &
           4*Mf*(2*Mf2*ml4 + (Mf4 + t2 + 2*s*(t - Mf2))*ml2 - &
           t*(2*s2 + 2*(t - Mf2)*s + (Mf2 + t)**2))*Mi + &
           ml4*(4*Mf4 - 4*t*Mf2) - ml2*(t - Mf2)*(3*Mf4 + t2 + &
           2*s*(t - 3*Mf2)) + t*(-3*Mf6 - t*Mf4 + 3*t2*Mf2 + t3 + &
           2*s2*(t - 3*Mf2) + 2*s*(3*Mf4 - 4*t*Mf2 + t2)))*c2))*C3a + &
           (1./(3*M4*Mf2))*(4*c2*C3v**2*(-2*M2*ml2*Mi6 - 2*M2*t*Mi6 - &
           4*M2*ml4*Mi4 + 6*M2*t2*Mi4 - 6*M2*ml2*Mf2*Mi4 - &
           6*M2*t*Mf2*Mi4 + 4*M2*ml2*s*Mi4 + 2*M2*ml2*t*Mi4 + &
           4*M2*s*t*Mi4 + 2*M2*ml2*Mf4*Mi2 - 2*M2*t*Mf4*Mi2 - &
           6*M2*t3*Mi2 + 2*M2*ml2*t2*Mi2 - 8*M2*s*t2*Mi2 - &
           8*M2*ml4*Mf2*Mi2 + 8*M2*t2*Mf2*Mi2 + 8*M2*ml2*s*Mf2*Mi2 + &
           16*M2*s*t*Mf2*Mi2 + 4*M2*ml4*t*Mi2 - 4*M2*s2*t*Mi2 - &
           8*M2*ml4*Mf3*Mi + 16*M2*t2*Mf3*Mi - 8*M2*ml2*t*Mf3*Mi + &
           6*M2*ml2*Mf6 - 6*M2*t*Mf6 + 2*M2*t4 + 4*M2*ml4*Mf4 + &
           6*M2*t2*Mf4 - 12*M2*ml2*s*Mf4 - 10*M2*ml2*t*Mf4 + 12*M2*s*t*Mf4 - &
           2*M2*ml2*t3 + 4*M2*s*t3 + 4*M2*s2*t2 - 4*M2*ml2*s*t2 - &
           2*M2*t3*Mf2 + 6*M2*ml2*t2*Mf2 - 16*M2*s*t2*Mf2 - &
           4*M2*ml4*t*Mf2 - 12*M2*s2*t*Mf2 + 16*M2*ml2*s*t*Mf2)) + &
           (1./(3*M4*Mf2))*(4*c2*C6a**2*(ml4*Mi6 - ml2*t*Mi6 + &
           2*ml4*Mf*Mi5 - 2*ml2*t*Mf*Mi5 + 3*ml2*t2*Mi4 - ml4*Mf2*Mi4 + &
           ml2*t*Mf2*Mi4 - 3*ml4*t*Mi4 - 4*ml4*Mf3*Mi3 + 4*ml2*t*Mf3*Mi3 + &
           4*ml2*t2*Mf*Mi3 - 4*ml4*t*Mf*Mi3 - ml4*Mf4*Mi2 + ml2*t*Mf4*Mi2 - &
           3*ml2*t3*Mi2 + 3*ml4*t2*Mi2 + 2*ml2*t2*Mf2*Mi2 - &
           2*ml4*t*Mf2*Mi2 + 2*ml4*Mf5*Mi - 2*ml2*t*Mf5*Mi + &
           4*ml2*t2*Mf3*Mi - 4*ml4*t*Mf3*Mi - 2*ml2*t3*Mf*Mi + &
           2*ml4*t2*Mf*Mi + ml4*Mf6 - ml2*t*Mf6 + ml2*t4 + 3*ml2*t2*Mf4 - &
           3*ml4*t*Mf4 - ml4*t3 - 3*ml2*t3*Mf2 + 3*ml4*t2*Mf2)) + &
           (1./(3*M4*Mf2))*(4*c2*C4a**2*(2*ml2*Mf8 - 2*t*Mf8 + 4*Mi*ml2*Mf7 - &
           4*Mi*t*Mf7 + 4*ml4*Mf6 + 2*Mi2*ml2*Mf6 - 2*t2*Mf6 - 4*ml2*s*Mf6 - &
           2*Mi2*t*Mf6 - 2*ml2*t*Mf6 + 4*s*t*Mf6 + 8*Mi*ml4*Mf5 - 8*Mi*t2*Mf5 - &
           8*Mi*ml2*s*Mf5 + 8*Mi*s*t*Mf5 + 2*t3*Mf4 - 2*Mi4*ml2*Mf4 + &
           2*ml2*t2*Mf4 - 8*s*t2*Mf4 - 2*Mi4*t*Mf4 - 4*ml4*t*Mf4 - &
           4*s2*t*Mf4 + 8*Mi2*s*t*Mf4 + 8*ml2*s*t*Mf4 - 8*Mi3*ml4*Mf3 - &
           4*Mi*t3*Mf3 - 4*Mi5*ml2*Mf3 + 8*Mi3*t2*Mf3 + 4*Mi*ml2*t2*Mf3 - &
           8*Mi*s*t2*Mf3 + 8*Mi3*ml2*s*Mf3 - 4*Mi5*t*Mf3 - 8*Mi*s2*t*Mf3 + &
           8*Mi3*s*t*Mf3 + 8*Mi*ml2*s*t*Mf3 - 4*Mi4*ml4*Mf2 + 2*t4*Mf2 - &
           6*Mi2*t3*Mf2 - 2*ml2*t3*Mf2 + 4*s*t3*Mf2 - 2*Mi6*ml2*Mf2 + &
           6*Mi4*t2*Mf2 + 2*Mi2*ml2*t2*Mf2 + 4*s2*t2*Mf2 - &
           8*Mi2*s*t2*Mf2 - 4*ml2*s*t2*Mf2 + 4*Mi4*ml2*s*Mf2 - 2*Mi6*t*Mf2 + &
           4*Mi2*ml4*t*Mf2 + 2*Mi4*ml2*t*Mf2 - 4*Mi2*s2*t*Mf2 + &
           4*Mi4*s*t*Mf2)) + (1./(3*M4*Mf2))*(4*c2*C4v**2* &
           (2*ml2*Mf8 - 2*t*Mf8 - 4*Mi*ml2*Mf7 + 4*Mi*t*Mf7 + 4*ml4*Mf6 + &
           2*Mi2*ml2*Mf6 - 2*t2*Mf6 - 4*ml2*s*Mf6 - 2*Mi2*t*Mf6 - &
           2*ml2*t*Mf6 + 4*s*t*Mf6 - 8*Mi*ml4*Mf5 + 8*Mi*t2*Mf5 + &
           8*Mi*ml2*s*Mf5 - 8*Mi*s*t*Mf5 + 2*t3*Mf4 - 2*Mi4*ml2*Mf4 + &
           2*ml2*t2*Mf4 - 8*s*t2*Mf4 - 2*Mi4*t*Mf4 - 4*ml4*t*Mf4 - &
           4*s2*t*Mf4 + 8*Mi2*s*t*Mf4 + 8*ml2*s*t*Mf4 + 8*Mi3*ml4*Mf3 + &
           4*Mi*t3*Mf3 + 4*Mi5*ml2*Mf3 - 8*Mi3*t2*Mf3 - 4*Mi*ml2*t2*Mf3 + &
           8*Mi*s*t2*Mf3 - 8*Mi3*ml2*s*Mf3 + 4*Mi5*t*Mf3 + 8*Mi*s2*t*Mf3 - &
           8*Mi3*s*t*Mf3 - 8*Mi*ml2*s*t*Mf3 - 4*Mi4*ml4*Mf2 + 2*t4*Mf2 - &
           6*Mi2*t3*Mf2 - 2*ml2*t3*Mf2 + 4*s*t3*Mf2 - 2*Mi6*ml2*Mf2 + &
           6*Mi4*t2*Mf2 + 2*Mi2*ml2*t2*Mf2 + 4*s2*t2*Mf2 - &
           8*Mi2*s*t2*Mf2 - 4*ml2*s*t2*Mf2 + 4*Mi4*ml2*s*Mf2 - 2*Mi6*t*Mf2 + &
           4*Mi2*ml4*t*Mf2 + 2*Mi4*ml2*t*Mf2 - 4*Mi2*s2*t*Mf2 + &
           4*Mi4*s*t*Mf2))



      ampl1=(1./(3*M4*Mf2))*(4*c2*C5v**2* &
           (2*ml2*Mf8 - 2*t*Mf8 - 4*Mi*ml2*Mf7 + 4*Mi*t*Mf7 + ml4*Mf6 + &
           2*Mi2*ml2*Mf6 + 6*t2*Mf6 - 4*ml2*s*Mf6 - 2*Mi2*t*Mf6 - &
           7*ml2*t*Mf6 + 4*s*t*Mf6 - 2*Mi*ml4*Mf5 - 8*Mi*t2*Mf5 + &
           8*Mi*ml2*s*Mf5 + 10*Mi*ml2*t*Mf5 - 8*Mi*s*t*Mf5 - Mi2*ml4*Mf4 - &
           6*t3*Mf4 - 2*Mi4*ml2*Mf4 + 4*Mi2*t2*Mf4 + 9*ml2*t2*Mf4 - &
           12*s*t2*Mf4 - 2*Mi4*t*Mf4 - 3*ml4*t*Mf4 - 3*Mi2*ml2*t*Mf4 - &
           4*s2*t*Mf4 + 8*Mi2*s*t*Mf4 + 12*ml2*s*t*Mf4 + 4*Mi3*ml4*Mf3 + &
           4*Mi*t3*Mf3 + 4*Mi5*ml2*Mf3 - 8*Mi*ml2*t2*Mf3 + 16*Mi*s*t2*Mf3 - &
           8*Mi3*ml2*s*Mf3 + 4*Mi5*t*Mf3 + 4*Mi*ml4*t*Mf3 - 4*Mi3*ml2*t*Mf3 + &
           8*Mi*s2*t*Mf3 - 8*Mi3*s*t*Mf3 - 16*Mi*ml2*s*t*Mf3 - Mi4*ml4*Mf2 + &
           2*t4*Mf2 - 2*Mi2*t3*Mf2 - 5*ml2*t3*Mf2 + 12*s*t3*Mf2 - &
           2*Mi6*ml2*Mf2 + 2*Mi4*t2*Mf2 + 3*ml4*t2*Mf2 + 8*s2*t2*Mf2 - &
           16*Mi2*s*t2*Mf2 - 12*ml2*s*t2*Mf2 + 4*Mi4*ml2*s*Mf2 - &
           2*Mi6*t*Mf2 + 2*Mi2*ml4*t*Mf2 + 7*Mi4*ml2*t*Mf2 - &
           4*Mi2*s2*t*Mf2 + 4*Mi4*s*t*Mf2 - 2*Mi5*ml4*Mf + 2*Mi*ml2*t3*Mf - &
           8*Mi*s*t3*Mf - 2*Mi*ml4*t2*Mf + 4*Mi3*ml2*t2*Mf - 8*Mi*s2*t2*Mf + &
           8*Mi3*s*t2*Mf + 8*Mi*ml2*s*t2*Mf - 4*Mi3*ml4*t*Mf - 6*Mi5*ml2*t*Mf + &
           8*Mi3*ml2*s*t*Mf + Mi6*ml4 + ml2*t4 - 4*s*t4 - ml4*t3 + &
           Mi2*ml2*t3 - 4*s2*t3 + 8*Mi2*s*t3 + 4*ml2*s*t3 - &
           Mi2*ml4*t2 - 5*Mi4*ml2*t2 + 4*Mi2*s2*t2 - 4*Mi4*s*t2 + &
           Mi4*ml4*t + 3*Mi6*ml2*t - 4*Mi4*ml2*s*t)) + &
           (16*c2*C5a**2*(Mi2 + 2*Mf*Mi + Mf2 - t)* &
           (ml4 - (-2*Mf2 + 2*s + t)*ml2 + s2 - s*Mf2 - 2*t*Mf2 + s*t + &
           Mi2*(ml2 + Mf2 - s)))/(3*Mf2) - &
           (16*c2*C5a*C6a*ml2*(Mi2 + 2*Mf*Mi + Mf2 - t)* &
           (Mi4 + (ml2 - Mf2 - s - 2*t)*Mi2 + t2 + s*Mf2 - t*Mf2 + s*t + &
           ml2*(Mf2 - t)))/(3*M2*Mf2) &
           + C4v*((1./(3*M4*Mf2))*(4*C5v*(4*ml2*Mf8 - 4*t*Mf8 - &
           8*Mi*ml2*Mf7 + 8*Mi*t*Mf7 + 4*ml4*Mf6 + 4*Mi2*ml2*Mf6 + &
           4*t2*Mf6 - 8*ml2*s*Mf6 - 4*Mi2*t*Mf6 - 8*ml2*t*Mf6 + 8*s*t*Mf6 - &
           8*Mi*ml4*Mf5 + 16*Mi*ml2*s*Mf5 + 8*Mi*ml2*t*Mf5 - 16*Mi*s*t*Mf5 + &
           4*t3*Mf4 - 4*Mi4*ml2*Mf4 + 4*ml2*t2*Mf4 - 16*s*t2*Mf4 - &
           4*Mi4*t*Mf4 - 8*ml4*t*Mf4 - 8*s2*t*Mf4 + 16*Mi2*s*t*Mf4 + &
           16*ml2*s*t*Mf4 + 8*Mi3*ml4*Mf3 - 8*Mi*t3*Mf3 + 8*Mi5*ml2*Mf3 + &
           16*Mi*s*t2*Mf3 - 16*Mi3*ml2*s*Mf3 + 8*Mi5*t*Mf3 + &
           8*Mi*ml4*t*Mf3 - 8*Mi3*ml2*t*Mf3 + 16*Mi*s2*t*Mf3 - &
           16*Mi3*s*t*Mf3 - 16*Mi*ml2*s*t*Mf3 - 4*Mi4*ml4*Mf2 - 4*t4*Mf2 + &
           4*Mi2*t3*Mf2 + 8*s*t3*Mf2 - 4*Mi6*ml2*Mf2 + 4*Mi4*t2*Mf2 + &
           4*ml4*t2*Mf2 - 4*Mi2*ml2*t2*Mf2 + 8*s2*t2*Mf2 - &
           16*Mi2*s*t2*Mf2 - 8*ml2*s*t2*Mf2 + 8*Mi4*ml2*s*Mf2 - &
           4*Mi6*t*Mf2 + 8*Mi4*ml2*t*Mf2 - 8*Mi2*s2*t*Mf2 + 8*Mi4*s*t*Mf2)* &
           c2) + (1./(3*M4*Mf2))*(4*C4a*(2*a*ml2*Mf8 - 2*a*t*Mf8 - &
           6*a*Mi2*ml2*Mf6 - 2*a*t2*Mf6 + 2*a*Mi2*t*Mf6 + 2*a*ml2*t*Mf6 + &
           4*a*s*t*Mf6 + 2*a*t3*Mf4 + 6*a*Mi4*ml2*Mf4 - 4*a*Mi2*t2*Mf4 - &
           2*a*ml2*t2*Mf4 + 8*a*s*t2*Mf4 + 2*a*Mi4*t*Mf4 - &
           4*a*Mi2*ml2*t*Mf4 - 8*a*Mi2*s*t*Mf4 + 2*a*t4*Mf2 - &
           6*a*Mi2*t3*Mf2 - 2*a*ml2*t3*Mf2 + 4*a*s*t3*Mf2 - &
           2*a*Mi6*ml2*Mf2 + 6*a*Mi4*t2*Mf2 + 2*a*Mi2*ml2*t2*Mf2 - &
           8*a*Mi2*s*t2*Mf2 - 2*a*Mi6*t*Mf2 + 2*a*Mi4*ml2*t*Mf2 + &
           4*a*Mi4*s*t*Mf2)*c2) + (16*a*C5a*(Mi2 - Mf2 - t)* &
           ((ml2 + t)*Mi2 + ml2*(t - Mf2) + t*(Mf2 - 2*s - t))*c2)/(3*M2) + &
           (8*a*C3a*(Mi2 - Mf2 - t)*(Mi2 - 4*Mf*Mi + 3*Mf2 - t)* &
           ((ml2 + t)*Mi2 + ml2*(t - Mf2) + t*(Mf2 - 2*s - t))*c2)/ &
           (3*M**3*Mf)) + C5v* &
           ((1./(3*M4*Mf2))*(4*C4a*(2*a*ml2*Mf8 - 2*a*t*Mf8 - &
           6*a*Mi2*ml2*Mf6 + 2*a*t2*Mf6 + 2*a*Mi2*t*Mf6 - 2*a*ml2*t*Mf6 + &
           4*a*s*t*Mf6 + 2*a*t3*Mf4 + 6*a*Mi4*ml2*Mf4 - 4*a*Mi2*t2*Mf4 - &
           2*a*ml2*t2*Mf4 + 2*a*Mi4*t*Mf4 + 4*a*Mi2*ml2*t*Mf4 - &
           8*a*Mi2*s*t*Mf4 - 2*a*t4*Mf2 + 2*a*Mi2*t3*Mf2 + 2*a*ml2*t3*Mf2 - &
           4*a*s*t3*Mf2 - 2*a*Mi6*ml2*Mf2 + 2*a*Mi4*t2*Mf2 + &
           2*a*Mi2*ml2*t2*Mf2 - 2*a*Mi6*t*Mf2 - 2*a*Mi4*ml2*t*Mf2 + &
           4*a*Mi4*s*t*Mf2)*c2) + (16*a*C5a*(Mi2 - Mf2 + t)* &
           ((ml2 + t)*Mi2 + ml2*(t - Mf2) + t*(Mf2 - 2*s - t))*c2)/(3*M2) + &
           (8*a*C3a*(Mi2 - Mf2 + t)*(Mi2 - 4*Mf*Mi + 3*Mf2 - t)* &
           ((ml2 + t)*Mi2 + ml2*(t - Mf2) + t*(Mf2 - 2*s - t))*c2)/ &
           (3*M**3*Mf))


      ampl2= C3v* &
           ((1./(3*M4*Mf2))*(4*C4v*(6*M*ml2*Mf7 - 6*M*t*Mf7 - &
           8*M*Mi*ml2*Mf6 + 8*M*Mi*t*Mf6 + 8*M*ml4*Mf5 + 2*M*Mi2*ml2*Mf5 - &
           2*M*t2*Mf5 - 12*M*ml2*s*Mf5 - 2*M*Mi2*t*Mf5 - 6*M*ml2*t*Mf5 + &
           12*M*s*t*Mf5 - 16*M*Mi*ml4*Mf4 + 16*M*Mi*t2*Mf4 + &
           16*M*Mi*ml2*s*Mf4 - 16*M*Mi*s*t*Mf4 - 4*M*Mi2*ml4*Mf3 + &
           6*M*t3*Mf3 - 6*M*Mi4*ml2*Mf3 + 2*M*ml2*t2*Mf3 - 16*M*s*t2*Mf3 + &
           8*M*Mi2*ml2*s*Mf3 - 6*M*Mi4*t*Mf3 - 8*M*ml4*t*Mf3 + &
           4*M*Mi2*ml2*t*Mf3 - 12*M*s2*t*Mf3 + 16*M*Mi2*s*t*Mf3 + &
           16*M*ml2*s*t*Mf3 + 16*M*Mi3*ml4*Mf2 + 8*M*Mi*t3*Mf2 + &
           8*M*Mi5*ml2*Mf2 - 16*M*Mi3*t2*Mf2 - 8*M*Mi*ml2*t2*Mf2 + &
           16*M*Mi*s*t2*Mf2 - 16*M*Mi3*ml2*s*Mf2 + 8*M*Mi5*t*Mf2 + &
           16*M*Mi*s2*t*Mf2 - 16*M*Mi3*s*t*Mf2 - 16*M*Mi*ml2*s*t*Mf2 - &
           4*M*Mi4*ml4*Mf + 2*M*t4*Mf - 6*M*Mi2*t3*Mf - 2*M*ml2*t3*Mf + &
           4*M*s*t3*Mf - 2*M*Mi6*ml2*Mf + 6*M*Mi4*t2*Mf + 2*M*Mi2*ml2*t2*Mf + &
           4*M*s2*t2*Mf - 8*M*Mi2*s*t2*Mf - 4*M*ml2*s*t2*Mf + &
           4*M*Mi4*ml2*s*Mf - 2*M*Mi6*t*Mf + 4*M*Mi2*ml4*t*Mf + &
           2*M*Mi4*ml2*t*Mf - 4*M*Mi2*s2*t*Mf + 4*M*Mi4*s*t*Mf)*c2) + &
           (1./(3*M4*Mf2))*(4*C5v*(6*M*ml2*Mf7 - 6*M*t*Mf7 - &
           8*M*Mi*ml2*Mf6 + 8*M*Mi*t*Mf6 + 4*M*ml4*Mf5 + 2*M*Mi2*ml2*Mf5 + &
           10*M*t2*Mf5 - 12*M*ml2*s*Mf5 - 2*M*Mi2*t*Mf5 - 14*M*ml2*t*Mf5 + &
           12*M*s*t*Mf5 - 8*M*Mi*ml4*Mf4 + 16*M*Mi*ml2*s*Mf4 + &
           8*M*Mi*ml2*t*Mf4 - 16*M*Mi*s*t*Mf4 - 4*M*Mi2*ml4*Mf3 - &
           2*M*t3*Mf3 - 6*M*Mi4*ml2*Mf3 + 10*M*ml2*t2*Mf3 - 24*M*s*t2*Mf3 + &
           8*M*Mi2*ml2*s*Mf3 - 6*M*Mi4*t*Mf3 - 8*M*ml4*t*Mf3 + &
           4*M*Mi2*ml2*t*Mf3 - 12*M*s2*t*Mf3 + 16*M*Mi2*s*t*Mf3 + &
           24*M*ml2*s*t*Mf3 + 8*M*Mi3*ml4*Mf2 - 8*M*Mi*t3*Mf2 + &
           8*M*Mi5*ml2*Mf2 + 16*M*Mi*s*t2*Mf2 - 16*M*Mi3*ml2*s*Mf2 + &
           8*M*Mi5*t*Mf2 + 8*M*Mi*ml4*t*Mf2 - 8*M*Mi3*ml2*t*Mf2 + &
           16*M*Mi*s2*t*Mf2 - 16*M*Mi3*s*t*Mf2 - 16*M*Mi*ml2*s*t*Mf2 - &
           2*M*t4*Mf + 2*M*Mi2*t3*Mf - 2*M*ml2*t3*Mf + 12*M*s*t3*Mf - &
           2*M*Mi6*ml2*Mf + 2*M*Mi4*t2*Mf + 4*M*ml4*t2*Mf - &
           6*M*Mi2*ml2*t2*Mf + 12*M*s2*t2*Mf - 16*M*Mi2*s*t2*Mf - &
           12*M*ml2*s*t2*Mf + 4*M*Mi4*ml2*s*Mf - 2*M*Mi6*t*Mf + &
           4*M*Mi2*ml4*t*Mf + 10*M*Mi4*ml2*t*Mf - 4*M*Mi2*s2*t*Mf + &
           4*M*Mi4*s*t*Mf - 8*M*Mi2*ml2*s*t*Mf)*c2) + &
           (1./(3*M4*Mf2))*(4*C4a*(6*a*M*ml2*Mf7 - 6*a*M*t*Mf7 + &
           8*a*M*Mi*ml2*Mf6 - 8*a*M*Mi*t*Mf6 - 10*a*M*Mi2*ml2*Mf5 + &
           2*a*M*t2*Mf5 - 2*a*M*Mi2*t*Mf5 - 2*a*M*ml2*t*Mf5 + 12*a*M*s*t*Mf5 - &
           16*a*M*Mi3*ml2*Mf4 + 16*a*M*Mi*s*t*Mf4 + 6*a*M*t3*Mf3 + &
           2*a*M*Mi4*ml2*Mf3 - 12*a*M*Mi2*t2*Mf3 - 6*a*M*ml2*t2*Mf3 + &
           8*a*M*s*t2*Mf3 + 6*a*M*Mi4*t*Mf3 + 4*a*M*Mi2*ml2*t*Mf3 - &
           8*a*M*Mi2*s*t*Mf3 + 8*a*M*Mi*t3*Mf2 + 8*a*M*Mi5*ml2*Mf2 - &
           16*a*M*Mi3*t2*Mf2 - 8*a*M*Mi*ml2*t2*Mf2 + 16*a*M*Mi*s*t2*Mf2 + &
           8*a*M*Mi5*t*Mf2 - 16*a*M*Mi3*s*t*Mf2 - 2*a*M*t4*Mf + &
           6*a*M*Mi2*t3*Mf + 2*a*M*ml2*t3*Mf - 4*a*M*s*t3*Mf + &
           2*a*M*Mi6*ml2*Mf - 6*a*M*Mi4*t2*Mf - 2*a*M*Mi2*ml2*t2*Mf + &
           8*a*M*Mi2*s*t2*Mf + 2*a*M*Mi6*t*Mf - 2*a*M*Mi4*ml2*t*Mf - &
           4*a*M*Mi4*s*t*Mf)*c2) - (16*a*C5a*(Mi2 + 4*Mf*Mi + 3*Mf2 - t)* &
           ((ml2 + t)*Mi2 + ml2*(t - Mf2) + t*(Mf2 - 2*s - t))*c2)/(3*M*Mf) + &
           (16*a*C3a*(Mi4 - 2*(t - Mf2)*Mi2 - 3*Mf4 + t2)* &
           ((ml2 + t)*Mi2 + ml2*(t - Mf2) + t*(Mf2 - 2*s - t))*c2)/ &
           (3*M2*Mf2)) + C4a* &
           ((1./(3*M4*Mf2))*(4*C6a*(-4*ml4*Mf6 + 4*ml2*t*Mf6 - &
           8*Mi*ml4*Mf5 + 8*Mi*ml2*t*Mf5 - 8*ml2*t2*Mf4 + 8*ml4*t*Mf4 + &
           8*Mi2*ml2*t*Mf4 - 8*ml2*s*t*Mf4 + 8*Mi3*ml4*Mf3 - &
           8*Mi*ml2*t2*Mf3 + 8*Mi*ml4*t*Mf3 + 8*Mi3*ml2*t*Mf3 - &
           16*Mi*ml2*s*t*Mf3 + 4*Mi4*ml4*Mf2 + 4*ml2*t3*Mf2 - &
           4*ml4*t2*Mf2 - 8*Mi2*ml2*t2*Mf2 + 8*ml2*s*t2*Mf2 + &
           4*Mi4*ml2*t*Mf2 - 8*Mi2*ml2*s*t*Mf2)*c2) + &
           (32*C5a*(Mi2 + 2*Mf*Mi + Mf2 - t)*(ml4 + (Mf2 - s)*ml2 - &
           t*(-Mi2 + Mf2 + t))*c2)/(3*M2))


      !differential cross section

      nu_MaElDELTA=(ampl0+ampl1+ampl2)*coupling


    end function nu_MaElDELTA


  end function nuMaEl


end module NeutrinoMatrixElement
