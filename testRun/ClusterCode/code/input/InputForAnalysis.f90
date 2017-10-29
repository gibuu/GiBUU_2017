!***************************************************************************
!****m* /inputForAnalysis
! NAME
! module inputForAnalysis
! PURPOSE
! This is the main module for reading some parameters concerning 
! the analysis routines
!***************************************************************************
module inputForAnalysis

  !*************************************************************************
  !****g* inputForAnalysis/yproj
  ! SOURCE
  !
  Real,save :: yproj=1.
  !
  ! PURPOSE
  ! rapidity of projectile in the CM frame of the two nucleus.
  ! NOTES
  ! Important only for heavy-ion collisions.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputForAnalysis/u0_Proj
  ! SOURCE
  !
  Real,save :: u0_Proj=1.
  !
  ! PURPOSE
  ! gamma factor of projectile in the CM frame of the two nucleus.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputForAnalysis/uz_Proj
  ! SOURCE
  !
  Real,save :: uz_Proj=1.
  !
  ! PURPOSE
  ! gamma*beta factor of projectile in the CM frame of the two nucleus.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputForAnalysis/SpectCut
  ! SOURCE
  !
  Real,save :: SpectCut=1.
  !
  ! PURPOSE
  ! Cut for selecting only spectator matter.
  ! NOTES 
  ! Not freqently used. Important Only for Heavy-ion collisions.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputForAnalysis/ImpactParameter_bin
  ! SOURCE
  !
  Real,save :: ImpactParameter_bin=1.
  !
  ! PURPOSE
  ! Gives the impact parameter bin, in units of [fm].
  ! NOTES 
  ! Only important for p+X reactions, where all observables are integrated 
  ! over the impact parameter from b=0 up to b_max 
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputForAnalysis/ThetaMax
  ! SOURCE
  !
  Integer,save :: ThetaMax=16
  !
  ! PURPOSE
  ! Number of considered polar angles for double differential spectra
  ! NOTES 
  ! Only important for p+X reactions.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputForAnalysis/Angles
  ! SOURCE
  !
  Real, dimension(1:16) ,save :: Angles=1000.
  !
  ! PURPOSE
  ! Polar angles (read from jobCard).
  ! NOTES 
  ! Only important for p+X reactions.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputForAnalysis/deltaTheta
  ! SOURCE
  !
  Real ,save :: deltaTheta=5.
  !
  ! PURPOSE
  ! Polar angle bin, i.e. [Angles-deltaTheta,Angles+deltaTheta]
  ! NOTES 
  ! Only important for p+X reactions.
  !
  !*************************************************************************

contains

  !*************************************************************************
  !****s* InputforAnalysis/GetAnalysisParameters
  ! NAME
  ! subroutine GetAnalysisParameters
  !
  ! PURPOSE
  ! * reads NAMELIST/InputAnalysis from JobCard:
  ! * boost parameters, projectile rapidity in cms and cut for 
  !   separating spectator from participant matter.
  !*************************************************************************
  subroutine GetAnalysisParameters
    use inputSMM, only : EventType
    implicit none
    integer :: ios
    !-----------------------------------------------------------------------
    NAMELIST /InputAnalysis/u0_Proj,uz_Proj,SpectCut,ImpactParameter_bin
    rewind(5)
    read(5,nml=inputAnalysis,IOSTAT=ios)
    !-----------------------------------------------------------------------
    if (EventType==1) then !HIC
       yProj = 0.5*log( (u0_Proj+uz_Proj)/(u0_Proj-uz_Proj) )
    else
       yProj = 1.0
    endif
    !-----------------------------------------------------------------------
    if (EventType==1) then
       write(*,999)
       write(*,1000) u0_Proj
       write(*,1001) uz_Proj
       write(*,1002) yProj
       write(*,1003) SpectCut
       write(*,1004) ImpactParameter_bin
       write(*,9999)
    else
       write(*,999)
       write(*,1003) SpectCut
       write(*,1004) ImpactParameter_bin
       write(*,9999)
    endif
    !-----------------------------------------------------------------------
999  format(//,78('='),/,20('='),' START LISTING OF ANALYSIS PARAMETERS ',20('='),/,78('='))
1000 format('  * Boost parameter u0_Proj              = ',f8.4)
1001 format('  * Boost parameter uz_Proj              = ',f8.4)
1002 format('  * Projectile rapidity in CMS           = ',f8.4)
1003 format('  * Cut for spectator matter             = ',f8.4)
1004 format('  * Impact parameter bin (only for p+Au) = ',f8.4)
9999 format(78('='),/,20('='),'  END LISTING OF ANALYSIS PARAMETERS  ',20('='),/,78('='),//)

!***************************************************************************
  end subroutine GetAnalysisParameters  !***********************************
!***************************************************************************

  !*************************************************************************
  !****s* InputforAnalysis/GetParametersSpectra
  ! NAME
  ! subroutine GetParametersSpectra
  !
  ! PURPOSE
  ! * reads NAMELIST/InputForSpectra from JobCard:
  ! * Parameters for double differential cross sections
  !*************************************************************************
  subroutine GetParametersSpectra
    implicit none
    integer :: ios
    !-----------------------------------------------------------------------
    NAMELIST /InputForSpectra/Angles, deltaTheta
    rewind(5)
    read(5,nml=InputForSpectra,IOSTAT=ios)

!    write(*,*) Angles


!***************************************************************************
  end subroutine GetParametersSpectra  !************************************
!***************************************************************************


end module inputForAnalysis
