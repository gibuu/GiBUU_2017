program CompareDIS

  use constants
  use inputGeneral, only: readinputGeneral
  use particleDefinition
  use particleProperties
  use eN_eventDefinition
  use eN_event

  use hadronFormation, only : forceInitFormation
  use CollTools, only: SetSomeDefaults_PY


  implicit none

  Type(electronNucleon_event), save :: eNev1, eNev2

  real :: eps1 = 0.05
  real :: eps2 = 0.99
  real :: W, Q2
  logical :: flagOK

  real :: E_,nu_,Q2_,W_
  real :: fT1, fT2
  real :: fak1, fak2
  
!  integer, parameter :: nPart=15
  integer, parameter :: nPart=45

  type(particle),dimension(1:nPart) :: OutPart_DIS

  real :: XS1, XS2
  integer, parameter :: nEv = 100
  integer :: iEv,iQ


  
  call readinputGeneral
!  call init_database
  call forceInitFormation
  call InitParticleProperties

  call SetSomeDefaults_PY

  W = 1.9

!  call FixedVals
  
  call TabulateElectron
!  call TabulateElectronEps
  call TabulateNeutrino
!  call ListEventNeutrino
!  call TabulateAnalytic
  call TabulateAnalyticPythia
  
contains

  subroutine FixedVals

    use eventGenerator_eN_lowEnergy, only: init_DIS


    Q2 = 2.0
    W = 100

    call eNev_setProcess(eNev1, 1,1)
    call eNev_init_eWQ(eNev1, eps1, W, Q2, flagOK)

    call init_DIS(eNev1,OutPart_DIS,XS1)    

    write(*,*) W, Q2, eps1, XS1
    
    
    call eNev_setProcess(eNev2, 1,1)    
    call eNev_init_eWQ(eNev2, eps2, W, Q2, flagOK)

    call init_DIS(eNev2,OutPart_DIS,XS2)

    write(*,*) W, Q2, eps2, XS2

    ! we have checked:
    ! ouput in 'fort.113' from "CollectXS_class" gives e.g.
    !     === xSect === (corrected for cuts)
    !     VMD * hadron                :  3.051D-03 mb     0.5961E-03
    !     direct * hadron             :  7.893D-03 mb     0.5961E-03
    !     anomalous * hadron          :  2.482D-03 mb     0.5961E-03
    !     DIS * hadron                :  9.126D-03 mb     0.5961E-03
    !     SUM                         :  2.255D-02 mb  
    !     === xSect === (corrected for cuts)
    !     VMD * hadron                :  4.944D-03 mb     0.7690E-02
    !     direct * hadron             :  9.617D-03 mb     0.7690E-02
    !     anomalous * hadron          :  3.916D-03 mb     0.7690E-02
    !     DIS * hadron                :  1.464D-02 mb     0.7690E-02
    !     SUM                         :  3.312D-02 mb  
    ! value 'DIS * hadron' for eps=~1 corresponds identically with values
    ! given in [Friberg:2000ra, fig.8.a]
    
    stop
    
  end subroutine FixedVals
  
  subroutine TabulateElectron

    use eventGenerator_eN_lowEnergy, only: init_DIS

    real :: sum1,sum2
    real :: c1,c2, cc1,cc2, cca,ccb, x1,x2

    call eNev_setProcess(eNev1, 1,1)
    call eNev_setProcess(eNev2, 1,1)
    do iQ=1,45
!    do iQ=10,10
       Q2 = iQ * 0.1

       sum1 = 0
       sum2 = 0
       
       call eNev_init_eWQ(eNev1, eps1, W, Q2, flagOK)
       !     call write_electronNucleon_event(eNev1,.FALSE.,.TRUE.)
       E_ = eNev1%lepton_in%momentum(0)
       call eNeV_GetKinV(eNev1, nu_,Q2_,W_, fT=fT1)
       fak1 = 1e3*pi/(fT1*E_*(E_-nu_))
       x1 = Q2/(2*0.938*nu_)
       c1 = (1-x1)/(6.24)
       cca = (nu_/E_)**2 + 2*(1-Q2/(4*E_**2) - nu_/E_)/(1+Q2/nu_**2)
       ccb = (nu_/E_)**2 + 2*(1 - nu_/E_)
       cc1 = c1 * cca/ccb
       
       do iEv=1,nEv
          call init_DIS(eNev1,OutPart_DIS,XS1)
          sum1=sum1+XS1
       enddo

       call eNev_init_eWQ(eNev2, eps2, W, Q2, flagOK)
       !     call write_electronNucleon_event(eNev2,.FALSE.,.TRUE.)
       E_ = eNev2%lepton_in%momentum(0)
       call eNeV_GetKinV(eNev2, nu_,Q2_,W_, fT=fT2)
       fak2 = 1e3*pi/(fT2*E_*(E_-nu_))
       x2 = Q2/(2*0.938*nu_)
       c2 = (1-x2)/(6.24)
       cca = (nu_/E_)**2 + 2*(1-Q2/(4*E_**2) - nu_/E_)/(1+Q2/nu_**2)
       ccb = (nu_/E_)**2 + 2*(1 - nu_/E_)
       cc2 = c2 * cca/ccb
       
       do iEv=1,nEv
          call init_DIS(eNev2,OutPart_DIS,XS2)
          sum2=sum2+XS2
       enddo
       
!       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,XS1,XS2,fT1,fT2,fak1,fak2
       write(*,'(1P,20e13.5)') W,Q2,eps1,eps2,sum1/nEv,sum2/nEv,fT1,fT2,fak1,fak2,c1,c2,cc1,cc2,x1,x2
       write(110,'(1P,20e13.5)') W,Q2,eps1,eps2,sum1/nEv,sum2/nEv,fT1,fT2,fak1,fak2,c1,c2,cc1,cc2,x1,x2
       
    end do

  end subroutine TabulateElectron

  subroutine TabulateNeutrino

    use Coll_nuN, only: DoColl_nuN_Py

    real :: sum1,sum2
    real :: form1,form2

    call eNev_setProcess(eNev1, 1,1)
    call eNev_setProcess(eNev2, 1,1)
    do iQ=1,45
       Q2 = iQ * 0.1

       sum1 = 0
       sum2 = 0
       
       call eNev_init_eWQ(eNev1, eps1, W, Q2, flagOK)
!       call write_electronNucleon_event(eNev1,.FALSE.,.TRUE.)
       E_ = eNev1%lepton_in%momentum(0)
       call eNeV_GetKinV(eNev1, nu_,Q2_,W_, fT=fT1)
       fak1 = 1e3*pi/(fT1*E_*(E_-nu_))
!       form1 = 1.0
       form1 = Q2_**2/((Q2_+0.78**2)**2)
!       form1 = Q2_**2/((Q2_+0.78**2)**2)*(1-Q2_/(2*0.938*nu_))
       
       do iEv=1,nEv
!          call DoColl_nuN_Py(eNev1,outPart_DIS,flagOK,.false.,XS1)
          call DoColl_nuN_Py(eNev1,outPart_DIS,flagOK,.true.,XS1) ! massless
          sum1=sum1+XS1
       enddo
       
       call eNev_init_eWQ(eNev2, eps2, W, Q2, flagOK)
!       call write_electronNucleon_event(eNev2,.FALSE.,.TRUE.)
       E_ = eNev2%lepton_in%momentum(0)
       call eNeV_GetKinV(eNev2, nu_,Q2_,W_, fT=fT2)
       fak2 = 1e3*pi/(fT2*E_*(E_-nu_))
!       form2 = 1.0
       form2 = Q2_**2/((Q2_+0.78**2)**2)
!       form2 = Q2_**2/((Q2_+0.78**2)**2)*(1-Q2_/(2*0.938*nu_))
       
       do iEv=1,nEv
!          call DoColl_nuN_Py(eNev2,outPart_DIS,flagOK,.false.,XS2)
          call DoColl_nuN_Py(eNev2,outPart_DIS,flagOK,.true.,XS2) ! massless
          sum2=sum2+XS2
       enddo
       
!       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,XS1,XS2,fT1,fT2,fak1,fak2
!       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,sum1/nEv,sum2/nEv,fT1,fT2,fak1,fak2
       write(*,'(1P,20e13.5)') W,Q2,eps1,eps2,sum1/nEv,sum2/nEv,fT1,fT2,fak1,fak2,form1,form2
       write(111,'(1P,20e13.5)') W,Q2,eps1,eps2,sum1/nEv,sum2/nEv,fT1,fT2,fak1,fak2,form1,form2
       
    end do

  end subroutine TabulateNeutrino

  subroutine ListEventNeutrino

    use Coll_nuN, only: DoColl_nuN_Py

    Q2 = 1.0

    call eNev_setProcess(eNev2, 1,1)
    call eNev_init_eWQ(eNev2, eps2, W, Q2, flagOK)
!       call write_electronNucleon_event(eNev2,.FALSE.,.TRUE.)
    call eNeV_GetKinV(eNev2, nu_,Q2_,W_, fT=fT2)

       
!       do iEv=1,nEv
!          call DoColl_nuN_Py(eNev2,outPart_DIS,flagOK,.false.,XS2)
    call DoColl_nuN_Py(eNev2,outPart_DIS,flagOK,.true.,XS2)

!       enddo
       
  end subroutine ListEventNeutrino

  subroutine TabulateElectronEps

    use eventGenerator_eN_lowEnergy, only: init_DIS
    use minkowski
    
    real :: sum1,sum2

    integer :: iEps

    
    Q2 = 1.0

    do iEps=0,10

       if (iEps.eq.0) then
          eps1 = 0.05
       else if (iEps.eq.10) then
          eps1 = 0.99
       else
          eps1 = iEps*0.1
       endif
    
       sum1 = 0
       sum2 = 0
       fak2 = 0
       fT2 = 0
       
       call eNev_setProcess(eNev1, 1,1)
       call eNev_init_eWQ(eNev1, eps1, W, Q2, flagOK)
       call write_electronNucleon_event(eNev1,.FALSE.,.TRUE.)
       write(*,*) 's = ',abs4sq(eNev1%lepton_in%momentum+eNev1%nucleon_free%momentum)
       write(*,*) 'light-x = ',eNeV_Get_LightX(eNev1)

       
       E_ = eNev1%lepton_in%momentum(0)
       call eNeV_GetKinV(eNev1, nu_,Q2_,W_, fT=fT1)
       write(*,*) ' x      = ',Q2/(2*0.938*nu_)
       write(*,*) ' y      = ',nu_/E_
       
       fak1 = 1e3*pi/(fT1*E_*(E_-nu_))
       
       do iEv=1,nEv
          call init_DIS(eNev1,OutPart_DIS,XS1)
          sum1=sum1+XS1
       enddo
       
       call pystat(1)
       
!       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,XS1,XS2,fT1,fT2,fak1,fak2
       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,sum1/nEv,sum2/nEv,fT1,fT2,fak1,fak2
       write(123,'(1P,12e13.5)') W,Q2,eps1,eps2,sum1/nEv,sum2/nEv,fT1,fT2,fak1,fak2
       
    end do

  end subroutine TabulateElectronEps

  
  subroutine TabulateAnalytic

    use Coll_nuN, only: AnaEstimate

    real :: sum1,sum2
    real :: form1,form2

    call eNev_setProcess(eNev1, 1,1)
    call eNev_setProcess(eNev2, 1,1)
    do iQ=1,30
       Q2 = iQ * 0.1

       sum1 = 0
       sum2 = 0
       
       call eNev_init_eWQ(eNev1, eps1, W, Q2, flagOK)
!       call write_electronNucleon_event(eNev1,.FALSE.,.TRUE.)
       E_ = eNev1%lepton_in%momentum(0)
       call eNeV_GetKinV(eNev1, nu_,Q2_,W_, fT=fT1)
       fak1 = 1e3*pi/(fT1*E_*(E_-nu_))
!       fak1 = 1e3/(2*pi*fT1)
!       form1 = 1.
       form1 = Q2_**2/((Q2_+0.78**2)**2)
!       form1 = Q2_**2/((Q2_+0.78**2)**2)*(1-Q2_/(2*0.938*nu_))
       
       XS1 = AnaEstimate(eNev1) ! = dsigma/dE'dcost in mb/GeV
       
       call eNev_init_eWQ(eNev2, eps2, W, Q2, flagOK)
!       call write_electronNucleon_event(eNev2,.FALSE.,.TRUE.)
       E_ = eNev2%lepton_in%momentum(0)
       call eNeV_GetKinV(eNev2, nu_,Q2_,W_, fT=fT2)
       fak2 = 1e3*pi/(fT2*E_*(E_-nu_))
!       fak2 = 1e3/(2*pi*fT2)
!       form2 = 1.
       form2 = Q2_**2/((Q2_+0.78**2)**2)
!       form2 = Q2_**2/((Q2_+0.78**2)**2)*(1-Q2_/(2*0.938*nu_))
       
       XS2 = AnaEstimate(eNev2) ! = dsigma/dE'dcost in mb/GeV
       
       ! fak1*XS1 = 1/Gamma dsigma/dE'dOmega' in 
       
!       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,XS1,XS2,fT1,fT2,fak1,fak2
!       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,sum1/nEv,sum2/nEv,fT1,fT2,fak1,fak2
       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,form1*XS1,form2*XS2,fT1,fT2,fak1,fak2
       write(113,'(1P,12e13.5)') W,Q2,eps1,eps2,form1*XS1,form2*XS2,fT1,fT2,fak1,fak2
       write(123,'(1P,12e13.5)') W,Q2,eps1,eps2,XS1,XS2,fT1,fT2,fak1,fak2
    end do

  end subroutine TabulateAnalytic

  subroutine TabulateAnalyticPythia

    use Coll_nuN, only: AnaEstimatePythia

    real :: sum1,sum2
    real :: form1,form2

    call eNev_setProcess(eNev1, 1,1)
    call eNev_setProcess(eNev2, 1,1)
    do iQ=1,30
       Q2 = iQ * 0.1

       sum1 = 0
       sum2 = 0
       
       call eNev_init_eWQ(eNev1, eps1, W, Q2, flagOK)
!       call write_electronNucleon_event(eNev1,.FALSE.,.TRUE.)
       E_ = eNev1%lepton_in%momentum(0)
       call eNeV_GetKinV(eNev1, nu_,Q2_,W_, fT=fT1)
       fak1 = 1e3*pi/(fT1*E_*(E_-nu_))
!       fak1 = 1e3/(2*pi*fT1)
       form1 = 1.
!       form1 = 1./(1-Q2_/(2*0.938*nu_))
!       form1 = Q2_**2/((Q2_+0.78**2)**2)
!       form1 = Q2_**2/((Q2_+0.78**2)**2)*(1-Q2_/(2*0.938*nu_))
       
       XS1 = AnaEstimatePythia(eNev1) ! = dsigma/dE'dcost in mb/GeV
       
       call eNev_init_eWQ(eNev2, eps2, W, Q2, flagOK)
!       call write_electronNucleon_event(eNev2,.FALSE.,.TRUE.)
       E_ = eNev2%lepton_in%momentum(0)
       call eNeV_GetKinV(eNev2, nu_,Q2_,W_, fT=fT2)
       fak2 = 1e3*pi/(fT2*E_*(E_-nu_))
!       fak2 = 1e3/(2*pi*fT2)
       form2 = 1.
!       form2 = 1./(1-Q2_/(2*0.938*nu_))
!       form2 = Q2_**2/((Q2_+0.78**2)**2)
!       form2 = Q2_**2/((Q2_+0.78**2)**2)*(1-Q2_/(2*0.938*nu_))
       
       XS2 = AnaEstimatePythia(eNev2) ! = dsigma/dE'dcost in mb/GeV
       
       ! fak1*XS1 = 1/Gamma dsigma/dE'dOmega' in 
       
!       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,XS1,XS2,fT1,fT2,fak1,fak2
!       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,sum1/nEv,sum2/nEv,fT1,fT2,fak1,fak2
       write(*,'(1P,12e13.5)') W,Q2,eps1,eps2,form1*XS1,form2*XS2,fT1,fT2,fak1,fak2
       write(114,'(1P,12e13.5)') W,Q2,eps1,eps2,form1*XS1,form2*XS2,fT1,fT2,fak1,fak2
       
    end do

  end subroutine TabulateAnalyticPythia
    
end program CompareDIS
