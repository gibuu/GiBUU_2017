!***************************************************************************
!****m* /WriteStatus
! NAME
! module WriteStatus
!
! PURPOSE
! This module contains all routines involving information written on 
! the terminal.
!***************************************************************************
module WriteStatus
  
  PRIVATE 
  
  PUBLIC :: IOControl,WriteHypInfo
  
contains

  !*************************************************************************
  subroutine IOControl(loc,status,routine,variable)
  !*************************************************************************
    implicit none
    integer,      intent(in) :: status,loc
    character(*), intent(in) :: routine,variable
    if(status /= 0) then
       write(*,*) 'from subroutine ',routine,' : '
       if (loc==1) then !allocations
          write(*,*) 'allocation of variable: ',variable ,'NOT successfull'
       endif
       if (loc==2) then!deallocations
          write(*,*) 'deallocation of variable: ',variable ,'NOT successfull'
       endif
       write(*,*) '!!! TERMINATION OF PROGRAM NOW !!!'
       STOP
    endif
  !*************************************************************************
  end subroutine IOControl
  !*************************************************************************



  !*************************************************************************
  subroutine WriteHypInfo(i,j,SubEvents,NumEnsemples,Clusters, & 
       & Hyperons,Mesons,Lambdas,Sigmas,Pions,Kaonsp,Kaons0, & 
       & FragmentVector)
  !*************************************************************************
    use typeDefinitions, only : cluster
    implicit none
  !-------------------------------------------------------------------------
  ! Input Variables
  !-------------------------------------------------------------------------
    integer,                    intent(in) :: i,j,SubEvents,NumEnsemples
    integer,                    intent(in) :: Clusters,Hyperons,Mesons
    integer,                    intent(in) :: Lambdas,Sigmas,Kaonsp,Kaons0
    integer, dimension(-1:1),   intent(in) :: Pions
    type(cluster),dimension(:), intent(in) :: FragmentVector
  !-------------------------------------------------------------------------
  ! Local Variables
  !-------------------------------------------------------------------------
  integer       :: m
  integer, SAVE :: He3LYield,He3SYield,He3MYield
  integer, SAVE :: He4LYield,He4SYield,He4MYield
  integer, SAVE :: Hyp_Yield,Mes_Yield,HypC
  integer, SAVE :: Lambda_Yield,Sigma_Yield,Kaonp_Yield,Kaon0_Yield
  integer, dimension(-1:1), SAVE :: Pion_Yield

  !-------------------------------------------------------------------------
  if (i==1 .and. j==1) then
     He3LYield = 0
     He3SYield = 0
     He3MYield = 0
     He4LYield = 0
     He4SYield = 0
     He4MYield = 0
     Hyp_Yield = 0
     Mes_Yield = 0
     HypC      = 0
     Lambda_Yield = 0
     Sigma_Yield = 0
     Pion_Yield(-1:1) = 0
     Kaonp_Yield = 0
     Kaon0_Yield = 0
  endif
  !-------------------------------------------------------------------------
  do m=1,Clusters
     if (FragmentVector(m)%MassNumber.eq.3 .and. & 
          FragmentVector(m)%ChargeNumber.eq.2 ) then
        if (FragmentVector(m)%HypNumber.eq.1) then
           if (FragmentVector(m)%HypType=='L') then
              write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) '  >>>>>>>>>>>>>>> formation of He3Lambda! <<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) 'SubEvent, NumEnsemples, nvals,nucleons,Hyperons,Clusters  = ', & 
                   & i,j,m
              He3LYield = He3LYield + 1  !3HeL
           else
              write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) '  >>>>>>>>>>>>>>> formation of He3Sigma! <<<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) 'SubEvent, NumEnsemples, nvals,nucleons,Hyperons,Clusters  = ', & 
                   & i,j,m
              He3SYield = He3SYield + 1  !3HeL
           endif
        endif
        if ((FragmentVector(m)%HypNumber).eq.2) then
           He3MYield = He3MYield + 1
        endif
     endif

     if (FragmentVector(m)%MassNumber.eq.4 .and. & 
          FragmentVector(m)%ChargeNumber.eq.2 ) then
        if (FragmentVector(m)%HypNumber.eq.1) then
           if (FragmentVector(m)%HypType=='L') then
              write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) '  >>>>>>>>>>>>>>> formation of He4Lambda! <<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) 'SubEvent, NumEnsemples, nvals,nucleons,Hyperons,Clusters  = ', & 
                   & i,j,m
              He4LYield = He4LYield + 1  !3HeL
           else
              write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) '  >>>>>>>>>>>>>>> formation of He4Sigma! <<<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
              write(*,*) 'SubEvent, NumEnsemples, nvals,nucleons,Hyperons,Clusters  = ', & 
                   & i,j,m
              He4SYield = He4SYield + 1  !3HeL
           endif
        endif
        if ((FragmentVector(m)%HypNumber).eq.2) then
           He4MYield = He4MYield + 1
        endif
     endif
     HypC = HypC + FragmentVector(m)%HypNumber !Hyperons bounded in clusters
  end do

  Hyp_Yield        = Hyp_Yield        + Hyperons
  Mes_Yield        = Mes_Yield        + Mesons
  Lambda_Yield     = Lambda_Yield     + Lambdas
  Sigma_Yield      = Sigma_Yield      + Sigmas
  Pion_Yield(-1:1) = Pion_Yield(-1:1) + Pions(-1:1)
  Kaonp_Yield      = Kaonp_Yield      + Kaonsp
  Kaon0_Yield      = Kaon0_Yield      + Kaons0

  write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
  write(*,*) '  >>>>>>>>>>>>>>> formation of HyperClusters! <<<<<<<<<<<<<<<<<<<<< '
  write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '    
  write(*,*) 'SubEvent, NumEnsemples, nvals,nucleons,Hyperons,Clusters,FragmentVector  = ', & 
       i,j,clusters, & 
       '( ',(FragmentVector(m)%MassNumber, m=1,Clusters),' )', & 
       '( ',(FragmentVector(m)%HypNumber, m=1,Clusters),' )'
  !-------------------------------------------------------------------------
  if (i==SubEvents .and. j==NumEnsemples) then
     open(111,file='ParticleSpeciesResults.dat')
     write(111,997)
     write(111,998)
     write(111,999)
     write(111,1000) float(He3LYield)/float(SubEvents*NumEnsemples), & 
          & float(He3SYield)/float(SubEvents*NumEnsemples), & 
          & float(He3MYield)/float(SubEvents*NumEnsemples)
     write(111,1001) float(He4LYield)/float(SubEvents*NumEnsemples), & 
          & float(He4SYield)/float(SubEvents*NumEnsemples), & 
          & float(He4MYield)/float(SubEvents*NumEnsemples)
     write(111,997)
     write(111,1002)
     write(111,999)
     write(111,1003) float(Hyp_Yield) /float(SubEvents*NumEnsemples)
     write(111,1004) float(Lambda_Yield) / float(SubEvents*NumEnsemples), & 
          & float(Sigma_Yield) / float(SubEvents*NumEnsemples)
     write(111,1005) float(HypC) /float(SubEvents*NumEnsemples)
     write(111,1006) float(Mes_Yield) /float(SubEvents*NumEnsemples)
     write(111,1007) Pion_yield(-1:1) /float(SubEvents*NumEnsemples)
     write(111,1008) Kaon0_Yield/float(SubEvents*NumEnsemples), & 
          & Kaonp_Yield/float(SubEvents*NumEnsemples)
     close(111)
  endif
997  format('######################')
998  format('# Hypercluster Yields:')
999  format('######################')
1000 format('# He3L,He3S,He3M = ',3f15.6)
1001 format('# He4L,He4S,He4M = ',3f15.6)
1002 format('# Hypercluster Yields:')
1003 format('# Hyperons (L+S) = ',f15.6)
1004 format('# Lambdas,Sigmas = ',2f15.6)
1005 format('# Coala-Hyperons = ',f15.6)
1006 format('# Mesons (tot)   = ',f15.6)
1007 format('# \pi^{-,0,+}    = ',3f15.6)
1008 format('# K^{0,+}        = ',2f15.6)


  !*************************************************************************
  end subroutine WriteHypInfo
  !*************************************************************************


end module WriteStatus
