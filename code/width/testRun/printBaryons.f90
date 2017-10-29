program test
  use particleProperties
  use idTable
  use inputGeneral

  implicit none

  call init_Database
  call readinputGeneral





  call PrintBaryonProperties(lbound(baryon,dim=1),uBound(baryon,dim=1))

end program

  !*************************************************************************
  !****s* test/PrintBaryonProperties
  ! NAME
  ! subroutine PrintBaryonProperties(n,m,texFlag)
  ! PURPOSE
  ! given a starting index n and and an ending index m the field baryon(n:m)
  ! is printed to file "BaronProperties.tex"
  ! INPUTS
  ! * integer :: n,m               -- Array Indizes
  ! * logical, OPTIONAL :: texFlag -- to produce tex-Output in table form
  ! OUTPUT
  ! written to "BaryonProperties.tex","BaryonProperties_Decays.txt"
  !*************************************************************************
  subroutine PrintBaryonProperties(n,m,texFlag)
  use particleProperties
  use idTable
  use inputGeneral
    use IdTable
    use DecayChannels
    use baryonWidthVacuum

    implicit none
    integer, intent(in)          :: n,m
    logical, intent(in),optional :: texFlag

    integer :: j,zaehler,k
    character(75), parameter ::line='--------------------------------------------------------------------------'
    logical :: flag
    character(5), parameter :: stars='*****'
    logical :: once
 

    Open(100,file='BarDecays.tex')


    write(100,*) "% Replace _ by \_ to get a perfect result"
    write(100,*) "\begin{table} "
    write(100,*) "\begin{tabular} {l|lll} "
    write(100,'(" &   ",2(A20," & "),A25)') 'Particle \#1' , ' Particle \#2', 'Angular Momentum \\ '
    write(100,*) "\hline"
    Do k=lbound(baryon,dim=1),ubound(baryon,dim=1)
       once=.true.
       Do j=1,size(baryon(nucleon)%decays2Body)
          If(baryon(k)%decays2Body(j).gt.1E-06) then
             if(once) then
                Write(100,'("$",A20,"$&",F8.4,"&$ ",A12,"$ &$ ",A12,"$& ",I1,"\\")')baryon(k)%nameTex, baryon(k)%decays2Body(j) &
                     & ,meson(decays2Body_Baryon(j)%id(1))%nameTex &
                     & ,baryon(decays2Body_Baryon(j)%id(2))%nameTex &
                     & ,getAngularMomentum(j,k)
                once=.false.
             else
                Write(100,'("$",A20,"$&",F8.4,"&$ ",A12,"$ &$ ",A12,"$& ",I1,"\\")')"  ", baryon(k)%decays2Body(j) &
                     & ,meson(decays2Body_Baryon(j)%id(1))%nametex &
                     & ,baryon(decays2Body_Baryon(j)%id(2))%nametex &
                     & ,getAngularMomentum(j,k)
             end if
          end if
       end do
    End do
    write(100,*) "\end  tabular}"
    write(100,*) "\caption  2-body decay channels for the baryons in BUU}"
    write(100,*) "\end  table}"



  end subroutine PrintBaryonProperties

