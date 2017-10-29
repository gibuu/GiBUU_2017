program test
  use inputGeneral
  implicit none


  call readinputGeneral
  call init_database

!  call sub_photon
  call sub_electron

end program test




subroutine sub_electron
  ! Evaluates photon cross section by sigma(electron)/gamma
  use particleDefinition
  use idtable
  use degRad_conversion
  use resProd_lepton
  use vector, only : theta_in
  use constants, only : pi
  use output
  use photon_flux, only : gamma
  use minkowski, only : abs4,sp
  use particleProperties
  implicit none


  type(particle) :: p
  real :: baremass, xsection
  real, dimension(0:3) :: q,ein,eout
  real,parameter :: ebeam=1.5
  real,dimension(0:3) :: lout,pout
  real :: theta_elec
  real,parameter :: theta=32./180.*3.14
  real :: egamma,sigma
  integer :: resID,charge
  real :: xsec_spin_1_2,xsec_spin_5_2,xsec_spin_3_2,xsec
  real,dimension(delta:nres) :: xsec_per_res


  do charge=0,1
     p%charge=charge
     p%momentum=(/baryon(nucleon)%mass,0.,0.,0./)
     p%ID=nucleon

     if(charge.eq.1) then
        open(100,file='electron_p.dat')
     else
        open(100,file='electron_n.dat')
     end if
     write(100,*) '# E_beam=', Ebeam
     write(100,*) '# Theta [deg]=', degrees(Theta)
     write(100,'("# Column",I8,A40)') 1 , "Photon energy"
     write(100,'("# Column",I8,A40)') 2 , "Q^2"
     write(100,'("# Column",I8,A40)') 3 , "Total Cross section"
     write(100,'("# Column",I8,A40)') 4 , "Total Cross section, only spin 1/2"
     write(100,'("# Column",I8,A40)') 5 , "Total Cross section, only spin 3/2"
     write(100,'("# Column",I8,A40)') 6 , "Total Cross section, only spin 5/2"
     do resid=delta,nres
        write(100,'("# Column",I8,A40)') resid+5 , baryon(resID)%name
     end do

     do egamma=0.1,2,0.01
        if(egamma.gt.ebeam) exit
        sigma=0.

        xsec=0.
        xsec_spin_1_2=0.
        xsec_spin_3_2=0.
        xsec_spin_5_2=0.
        do resId=delta,nres
           eout=(/ebeam-egamma,0.,sin(theta)*(ebeam-egamma),cos(theta)*(ebeam-egamma)/)
           ein=(/ebeam,0.,0.,ebeam/)
           q=ein-eout
           theta_elec=theta_in(ein(1:3),eout(1:3))
           write(*,*) 'theta=', theta_elec
           xsec_per_res(resID)=dSigmadOmega_fdE_f_resProd(p,resID,theta_elec,ein(0),eout(0),q,lout,pout,bareMass)
           if(nint(baryon(resID)%spin*2.).eq.1) then
              xsec_spin_1_2=xsec_spin_1_2+xsec_per_res(resID)
           else if(nint(baryon(resID)%spin*2.).eq.3) then
              xsec_spin_3_2=xsec_spin_3_2+xsec_per_res(resID)
           else if(nint(baryon(resID)%spin*2.).eq.5) then
              xsec_spin_5_2=xsec_spin_5_2+xsec_per_res(resID)
           end if
        end do
        xsec=sum(xsec_per_res)

        write(100,'('//intTochar(2+nres-delta+5)//'E15.5)') egamma, -sp(q,q),xsec , xsec_spin_1_2,xsec_spin_3_2,xsec_spin_5_2,xsec_per_res


     end do
     close(100)
  end do
end subroutine sub_electron



subroutine sub_photon
  ! Evaluates directly photon induced resonance production cross section 
  ! Prints Xsections to files:
  ! * gamma_n_to_R_Xsection.dat
  ! * gamma_p_to_R_Xsection.dat
  ! * gamma_n_to_pipi_Xsection.dat
  ! * gamma_p_to_pipi_Xsection.dat
  ! * gamma_n_to_piN_Xsection.dat
  ! * gamma_p_to_piN_Xsection.dat

  use particleDefinition
  use idtable
  use resProd_lepton
  use output
  use particleProperties
  use electronPionProduction_medium
  use constants, only :pi
  implicit none

  type(particle) :: p
  real :: baremass_res, e_gamma,xsec
  real,dimension(delta:nres) :: xsec_per_res
  real,dimension(-1:1)::  xsection
  integer :: resid
  real, dimension(0:3) :: q,k,pf
  real :: xsec_spin_1_2,xsec_spin_5_2,xsec_spin_3_2,phi_k, theta_k, cosTheta
  real,parameter :: dCosTheta=0.01
  integer :: pionCharge,nucCharge
  logical :: success
  p%charge=0

  p%momentum=(/0.938,0.,0.,0./)
  p%mass=0.938
  p%ID=nucleon
  p%position=(/100.,100.,100./)

  ! Resonance contribution
  do nucCharge=0,1
     if(nucCharge.eq.0) then
        open(100,file='gamma_n_to_R_Xsection.dat')
        open(199,file='gamma_n_to_pipi_Xsection.dat')
     else
        open(100,file='gamma_p_to_R_Xsection.dat')
        open(199,file='gamma_p_to_pipi_Xsection.dat')
     end if
     p%charge=nucCharge


     write(100,'("# Column",I8,A40)') 1 , "Photon energy"
     write(100,'("# Column",I8,A40)') 2 , "Total Cross section"
     write(100,'("# Column",I8,A40)') 3 , "Total Cross section, only spin 1/2"
     write(100,'("# Column",I8,A40)') 4 , "Total Cross section, only spin 3/2"
     write(100,'("# Column",I8,A40)') 5 , "Total Cross section, only spin 5/2"

     do resid=delta,nres
        write(100,'("# Column",I8,A40)') resid+4 , baryon(resID)%name
     end do

     do e_gamma=0.1,1.4,0.01
        xsec=0.
        xsec_spin_1_2=0.
        xsec_spin_3_2=0.
        xsec_spin_5_2=0.
        do resid=delta,nres
           xsec_per_res(resID)=sigma_resProd(p,resId,(/e_gamma,0.,0.,e_gamma/),baremass_res)
           if(nint(baryon(resID)%spin*2.).eq.1) then
              xsec_spin_1_2=xsec_spin_1_2+xsec_per_res(resID)
           else if(nint(baryon(resID)%spin*2.).eq.3) then
              xsec_spin_3_2=xsec_spin_3_2+xsec_per_res(resID)
           else if(nint(baryon(resID)%spin*2.).eq.5) then
              xsec_spin_5_2=xsec_spin_5_2+xsec_per_res(resID)
           end if
        end do
        xsec=sum(xsec_per_res)

        write(100,'('//intTochar(2+nres-delta+4)//'E15.5)') e_gamma,  xsec , xsec_spin_1_2,xsec_spin_3_2,xsec_spin_5_2,xsec_per_res
        write(199,'(5E15.5)') e_gamma,  sigma_pipi_res(p,(/e_gamma,0.,0.,e_gamma/))
     end do

     close(100)
     close(199)
  end do


  do nucCharge=0,1
     if(nucCharge.eq.0) then
        open(66,file='gamma_n_to_piN_Xsection.dat')
     else
        open(66,file='gamma_p_to_piN_Xsection.dat')
     end if
     p%charge=nucCharge
     ! 1Pi contribution
     do e_gamma=0.,1.4,0.01
        xsection=0.
        write(*,*) 
        write(*,*) line
        write(*,*) 'E_gamma=',e_gamma
        do pionCharge=p%charge-1,p%charge
           q=(/e_gamma,0.,0.,e_gamma/)
           do cosTheta=-1,1,dCosTheta
              theta_k=acos(cosTheta)*180./pi
              phi_k=0.
              xsection(pionCharge)=xsection(pionCharge)+dSigmadOmega_k_med(p,pionCharge,phi_k,theta_k,q,k,pf,success) 
           end do
        end do
        xsection=2.*pi*dCosTheta*xsection
        write(66,'(5E15.4)') e_gamma,  xsection,sum(xsection)
     end do
     close(66)
  end do

end subroutine sub_photon
