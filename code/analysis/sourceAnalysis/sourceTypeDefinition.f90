!******************************************************************************
!****m* /sourceTypeDefinition
! NAME
! module sourceTypeDefinition
! PURPOSE
! Contains the type definition "quelle".
!******************************************************************************

module sourceTypeDefinition


  !****************************************************************************
  !****t* sourceTypeDefinition/Quelle
  ! NAME
  ! type Quelle
  ! PURPOSE
  ! Stores properties of fragmenting source into type "Quelle".
  !
  type quelle
     logical              :: status=.false.   !if true, source exists
     integer              :: Size=0     !size of existing source
     integer              :: Charge=0   !charge of existing source
     integer              :: nLambda=0  !number of Lambda-hyperons in the source
     integer              :: nSigma0=0  !number of Sigma0-hyperons in the source
     real, dimension(1:3) :: position=99999. !position of existing source [fm]
     real, dimension(1:3) :: velocity=99999. !average velocity of the source
     real                 :: ExEnergy=99999. !excitation energy of existing source [GeV/A]
     real                 :: radEnergy=99999. !energy of radial flow [GeV/A]
  end type quelle
  !
  !****************************************************************************

end module sourceTypeDefinition
