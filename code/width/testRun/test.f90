
program test
  use idTable
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: InitParticleProperties, hadron, nDecays
  implicit none

  write(*,*) "**********************************" 
  write(*,*) "**********************************" 
  write(*,*)
  Write(*,*)      "Testing the routines which generate the width of the resonances "
  write(*,*)
  write(*,*) "**********************************" 
  write(*,*) "**********************************" 
  write(*,*)

  call readinputGeneral
  call InitParticleProperties


 call testBaryon
!  call testDelta
!  call testDecay 
!   call testMeson


contains


subroutine testDelta
  use baryonWidth
  use baryonWidthMedium
  use mediumDefinition
  use constants, only : rhoNull

  type(medium) :: mediumAtPosition
  integer :: i, id,j
  real :: mass, mom(0:3)

  mom=(/0.,0.,0.25,0./)

  ! files : fort.100 = Full width, FullWidthMedium
  ! files : fort.200 = Sum of In-Width
  ! files : fort.300 = Sum of Out-Width
 
  mediumAtPosition%useMedium=.true.
  mediumAtPosition%densityProton=rhoNull/2.
  mediumAtPosition%densityNeutron=rhoNull/2
  mediumAtPosition%density=rhoNull


  write(*,*) "**********************************" 
  write(*,*)
  Write(*,*)      "Testing 'FullWidthBaryon for delta at rhoNull'"
  Open(999,file='DeltaWidth.dat')
  id=2
  
  do i=0,52  
     do j=0,20
        mom=(/0.,j*0.05,0.,0./)
        mass=1.071+float(i)*0.01
        write(999,'(6F15.9)') mass,mom(1), &
             & partialWidthBaryonMedium(id,mass,.false.,pion,nucleon,mom,mediumATposition)!, &
!             & mediumAtPosition%density, WidthBaryonMedium(id,mass,mom,  mediumAtPosition), &
!             & WidthBaryonMedium(id,mass,mom,  mediumAtPosition)- partialWidthBaryonMedium(id,mass,.false.,pion,nucleon,mom,mediumATposition)
     end do
  end do
  close(999)
end subroutine testDelta



!*************************************************************


  subroutine testDecay

  use baryonWidthMedium
  use mesonWidthMedium
  use mediumDefinition
!   use decayChannels

!   type(medium) :: mediumAtPosition
  integer :: i
  real :: mass
!   real :: mom(0:3)
  real, dimension(1:nDecays) :: decay  ! ,decayWidth
  logical ::  PauliFlag

!   mom=(/0.,0.,0.,0./)

  ! files : fort.100 = Full width, FullWidthMedium
  ! files : fort.200 = Sum of In-Width
  ! files : fort.300 = Sum of Out-Width

!   mediumAtPosition%useMedium=.false.
!   mediumAtPosition%densityProton=0.17

!   Allocate(decayWidth(1:size(decays2body_baryon)))

!!$
!!$  write(*,*) "**********************************" 
!!$  write(*,*)
!!$  Write(*,*)      "Testing Delta Decay Width'"
!!$  
!  mom=(/0.,0.,0.,0./)
!!$  mass=baryon(delta)%mass
!!$  call decayWidthBaryonMedium(delta,mass,mom,mediumATposition, decayWidth, pauliFlag)
!!$  Do i=1,size(decayWidth)
!!$     write(*,*) , i , decayWidth(i), pauliFlag
!!$  End Do

  mass=hadron(kaonStarBar)%mass-0.1
  Write(*,*)      "Testing Decay Width' of kaonStarBar"
  decay = decayWidthMesonMedium(kaonStarBar,mass,0,pauliFlag)
  Do i=1,nDecays
     write(*,*) i , decay(i), pauliFlag,mass
  End Do

  Write(*,*)      "Testing Decay Width' of kaonStar"
  mass=hadron(kaonStar)%mass-0.1
  decay = decayWidthMesonMedium(kaonStar,mass,0,pauliFlag)
  Do i=1,nDecays
     write(*,*) i , decay(i), pauliFlag,mass
  End Do
  
end subroutine testDecay



!***********************

  subroutine testBaryon
  use baryonWidth
  use baryonWidthMedium
  use mediumDefinition

  type(medium) :: mediumAtPosition
  integer :: i, id, mesonID, baryonID, index
  real :: mass,sum,feld(1:5),mom(0:3)

  mom=(/0.,0.,0.,0./)

  ! files : fort.100 = Full width, FullWidthMedium
  ! files : fort.200 = Sum of In-Width
  ! files : fort.300 = Sum of Out-Width


  write(*,*) "**********************************"
  write(*,*) "Testing 'FullWidthBaryon'"
  write(*,*) "**********************************"
  do id=1,nbar
     do i=1,400  
        mass = hadron(id)%minmass+i*0.01
        write(id+100,*) mass, FullWidthBaryon(id,mass), WidthBaryonMedium(id,mass,mom,mediumAtPosition), decayWidthBaryon(Id,mass)
     end do
  end do


  write(*,*) "******************************************"
  write(*,*) "Testing 'PartialWidthBaryon' for Out-Width"
  write(*,*) "******************************************"
  do id=1,nbar
     Print *, "resonance=",id
     do i=1,40  
        mass = hadron(id)%minmass+i*0.05
        sum=0.
        do mesonID=pion,dsStar_Minus
           do baryonID=nucleon,Omega_c
              sum=sum+partialwidthBaryon(ID,mass,.false.,mesonID,baryonID)
            end do
        end do
        If (Abs(sum-FullWidthBaryon(id,mass))>1E-5) then 
           Print *, "Sum of partial width not equal full width"
           Print *, "Resonance=",id
           Print *, "Mass=",mass
           Print *, sum, FullWidthBaryon(id,mass)
           stop
        end if
     end do
  end do
!!$

  write(*,*) "*******************************************"
  write(*,*) "Testing 'PartialWidthBaryon' for Out-Width"
  write(*,*) "Medium modified"
  write(*,*) "*******************************************"
  do id=1,nbar
     Print *, "resonance=",id
     do i=1,40  
        mass = hadron(id)%minmass+i*0.05
        sum=0.
        do mesonID=pion,dsStar_Minus
           do baryonID=nucleon,Omega_c
              sum=sum+partialwidthBaryonMedium(ID,mass,.false.,mesonID,baryonID,mom, mediumAtPosition)
            end do
        end do
        If (Abs(sum-widthBaryonMedium(id,mass,mom,MediumAtPosition))>1E-5) then 
           Print *, "Sum of partial width not equal full width"
           Print *, "Resonance=",id
           Print *, "Mass=",mass
           Print *, sum, FullWidthBaryon(id,mass), widthBaryonMedium(id,mass,mom,MediumAtPosition)
           stop
        end if
     end do
  end do



  write(*,*) "*****************************************"
  write(*,*) "Testing 'PartialWidthBaryon' for In-Width"
  write(*,*) "*****************************************"
  do id=1,nbar
     Print *, "resonance=",id
     write(200+id,*) "#Resonance ",id, "; Sum of In width over all baryon meson combinations"
     do i=1,50  
        mass = hadron(id)%minmass+i*0.04
        sum=0.
        feld=0.
        index=0
       do mesonID=pion,dsStar_Minus
           do baryonID=nucleon,Omega_c
              sum=sum+partialwidthBaryon(ID,mass,.true.,mesonID,baryonID)
              If (partialwidthBaryon(ID,mass,.true.,mesonID,baryonID)>0. .and. index<size(feld)) then
                 index=index+1
                 feld(index)=partialwidthBaryon(ID,mass,.true.,mesonID,baryonID)
              end if
           end do
        end do
        write(200+id,'(7F8.2)')  mass,sum,feld
      end do
   end do

   write(*,*) "******************************************"
   write(*,*) "Testing 'PartialWidthBaryon' for Out-Width"
   write(*,*) "******************************************"
   do id=1,nbar
      Print *, "resonance=",id
      write(300+id,*) "#Resonance ",id, "; Sum of In width over all baryon meson combinations"
      do i=1,50  
         mass = hadron(id)%minmass+i*0.04
         sum=0.
         feld=0.
         index=0
         do mesonID=pion,dsStar_Minus
            do baryonID=nucleon,Omega_c
               sum=sum+partialwidthBaryon(ID,mass,.false.,mesonID,baryonID)
               If (partialwidthBaryon(ID,mass,.false.,mesonID,baryonID)>0. .and. index<size(feld)) then
                  index=index+1
                  feld(index)=partialwidthBaryon(ID,mass,.false.,mesonID,baryonID)
               end if
            end do
         end do
         write(300+id,'(7F8.2)')  mass,sum,feld
      end do
   end do
 end subroutine testBaryon



 !****************************


 subroutine testMeson
   use mesonWidth
   use mesonWidthMedium
   use mediumDefinition

   integer :: i, id
   real :: mass,sum !,sum1
   type(medium) :: mediumAtPosition
   integer :: decayID_1,decayID_2,decayID_3
   integer :: charge1,charge2,charge3
   ! files : fort.5** = Full width
   ! files : fort.6** = Sum of partial Width
  real :: mom(0:3)

  mom=(/0.,0.,0.,0./)


   write(*,*) 
   write(*,*) "**********************************" 
   write(*,*) "**********************************" 
   write(*,*) 'TESTING MESON-SECTOR      '
   write(*,*) "**********************************" 
   write(*,*) "**********************************" 
   write(*,*) 


   write(*,*) "**********************************" 
   write(*,*)
   Write(*,*)      "Testing 'FullWidthMeson'"
   do id=pion,pion+nmes-1

      do i=1,10000  
         mass=i*0.0004
         write(id-pion+1+500,*) mass,FullWidthMeson(id,mass),WidthMesonMedium(id,mass,mom,mediumAtPosition)
      end do
   end do



   write(*,*) "**********************************" 
   write(*,*)
   Write(*,*)      "Testing 'Partial WidthMeson' "
   do id=pion,pion+nmes-1

      Print *, "resonance=",id
      write(id-pion+1+600,*) "#Resonance ",id, "; Sum of partial width "

      do i=0,100 
         mass=i*0.05
         sum=0.
         ! two Channels with only mesons
         do decayID_1=pion,dsStar_Minus
            do decayID_2=decayId_1,dsStar_Minus
               sum=sum+partialwidthMeson(ID,mass,decayID_1,decayID_2)
            end do
         end do
         ! two Channels with meson and photon
         do decayID_1=pion,dsStar_Minus
            sum=sum+partialwidthMeson(ID,mass,decayID_1,photon)
         end do
         sum=sum+partialwidthMeson(ID,mass,photon,photon)

!          sum1=sum
         ! three channels
         do decayID_1=pion,dsStar_Minus
            do decayID_2=decayID_1,dsStar_Minus
               do decayID_3=decayID_2,dsStar_Minus
                  do charge1=-1,1
                     do charge2=charge1,1
                        do charge3=charge2,1
!                            sum = sum + partialwidthMesonMedium(ID,mass,decayID_1,decayID_2,decayID_3,charge1,charge2,charge3, &
!                                                                mom,mediumAtPosition)
                           sum = sum + partialwidthMeson(ID,mass,decayID_1,decayID_2,decayID_3,charge1,charge2,charge3)
                        End do
                     End do
                  End do
               End do
            End do
         End do
         write(id-pion+1+600,'(3F12.8)')  mass,sum
      end do

   end do

  end subroutine testMeson



end program test
