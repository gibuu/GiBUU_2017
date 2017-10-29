!***************************************************************************
!****m* /typeDefinitions
! NAME
! module typeDefinitions
! FUNCTION
! This module contains all type definitions.
!***************************************************************************

module typeDefinitions

  !*************************************************************************
  !****t* typeDefinitions/Fragment
  ! NAME
  ! type Fragment
  ! PURPOSE
  ! This is the MAIN type definition for the fragment vector
  ! NOTES
  ! This type is common for both models of fragmentation (coalescence/SMM). 
  ! The element %momentumLRF is needed only in SMM-mode, in particular 
  ! for spectator fragmentation (ALADIN_Flag=true).
  !*************************************************************************
  Type cluster
     real, dimension (1:3)   :: position=0.    !position (fm)
     real, dimension (0:3)   :: momentum=0.    !4-momentum (GeV) in global CMS
     real, dimension (0:3)   :: momentumLRF=0. !4-momentum (GeV) in moving frame
     real                    :: mass=0.        !mass (GeV) (=A * M_part)
     integer                 :: ID=0           !ID (=1 for fragments)
     integer                 :: Mechanism=0    !Production mechanism (1:fission/2:evaporation)
     integer                 :: MassNumber=0   !mass number of fragment
     integer                 :: ChargeNumber=0 !charge number of fragment
     integer                 :: HypNumber=0    !strangeness content of fragment
     real                    :: pionic=0.      !Probability of hypernuclei from pion+Nuc
     character(2)            :: HypType='-'    !type of bounded hyperon(s)
     logical                 :: FreeBound=.true. !true: SMM-Cluster is nucleon, otherwise: fragment
     !stableFlag=.true. -> source decays into SMM clusters
     !stableFlag=.false.-> source's Ex < 0 ==> no decay via SMM. 
     logical                 :: stableFlag=.false. 
     integer, dimension(1:3) :: HypContent=-999     !index of bound hyperons
  End Type Cluster

  !*************************************************************************
  !****t* typeDefinitions/particle
  ! NAME
  ! type particle
  ! PURPOSE
  ! This is the MAIN type definition for the GiBUU particle vector
  ! NOTES
  ! Only real particles has to be considered!!!
  !*************************************************************************
  Type particle
     integer               :: number=0
     real                  :: bornTime=-999.
     real                  :: lastCollTime=-999.
     integer               :: collHis=-999
     real, dimension (1:3) :: position=0.
     real, dimension (0:3) :: momentum=0.
     real                  :: mass=0.
     integer               :: ID=999
     integer               :: Charge =0
     integer               :: event=0
     integer               :: ensemple=0
  End Type particle

  !*************************************************************************
  !****t* typeDefinitions/Quelle
  ! NAME
  ! type Quelle
  ! PURPOSE
  ! Stores properties of fragmenting source into type "Quelle".
  !
  type quelle
     logical              :: status=.false.   !if true, source exists
     integer              :: Size=0           !size of existing source
     integer              :: Charge=0         !charge of existing source
     real, dimension(1:3) :: position=99999.  !position of existing source [fm]
     real, dimension(1:3) :: velocity=99999.  !average velocity of the source
     real                 :: ExEnergy=99999.  !excitation energy of existing source [GeV/A]
     real                 :: radEnergy=99999. !energy of radial flow [GeV/A]
     integer              :: Type=99999       !defines source's origin (spectator/fireball)
  end type quelle
  !
  !*************************************************************************

end module typeDefinitions
