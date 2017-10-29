!******************************************************************************
!****m* /povray
! NAME
! module povray
! PURPOSE
! Includes routines to make povray output.
! Prerequisites:
! * POV-Ray (http://www.povray.org).
! * ImageMagick (http://www.imagemagick.org)
!******************************************************************************
module povray

  implicit none
  private

  public :: povray_output

contains

  !****************************************************************************
  !****s* povray/povray_output
  ! NAME
  ! subroutine povray_output (timestep, realPart, pertPart)
  ! INPUTS
  ! * integer :: timestep                                    --- current timestep number
  ! * type(particle), dimension(:,:) :: realPart, pertPart   --- real and perturbative particle vectors
  ! OUTPUT
  ! * at each timestep, three output files are written: Movie_Real_*.pov, Movie_Pert_*.pov, Movie_All_*.pov.
  ! PURPOSE
  ! Generates povray output files. Also a script "povray.sh" is generated,
  ! which needs to be executed when GIBUU is finished. This will then generate
  ! animated gif files out of the povray output.
  !****************************************************************************
  subroutine povray_output (timestep, realPart, pertPart)
    use particleDefinition
    use output, only: intTochar

    integer, intent(in) :: timestep
    type(particle), intent(in), dimension(:,:) :: realPart, pertPart

    call makePicture ('Movie_Real_'//intTochar(timestep)//'.pov', realPart)
    call makePicture ('Movie_Pert_'//intTochar(timestep)//'.pov', pertPart)
    call makePicture ('Movie_All_'//intTochar(timestep)//'.pov', realPart, pertPart)
  end subroutine


  !****************************************************************************
  !****s* povray/makePicture
  ! NAME
  ! subroutine makePicture (filename, particles, particles1)
  ! INPUTS
  ! * character(*)                             :: filename
  ! * type(particle), dimension(:,:)           :: particles
  ! * type(particle), dimension(:,:), optional :: particles1
  ! PURPOSE
  ! In file 'filename' a povray output according to the input vectors
  ! 'particles' and 'particles1' is generated. You can either give one or two
  ! particle vectors as input, since 'particles1' is optional.
  ! For each particle a sphere is generated around its position.
  ! The general settings and drawing macros are written to a
  ! povray include file ("povray_header.inc"). This can be modified after the
  ! run to adjust the viewing angle and drawing style etc.
  !****************************************************************************
  subroutine makePicture (filename, particles, particles1)
    use particleDefinition

    character(*), intent(in)                             :: filename
    type(particle), intent(in), dimension(:,:)           :: particles
    type(particle), intent(in), dimension(:,:), optional :: particles1

    integer :: i,j, nWritten
    logical, save :: initflag = .true.

    if (initFlag) then
      call writeHeader
      call writeScript
      initflag = .false.
    end if

    nWritten = 0
    open(77,file=filename)
    write(77,'(A)') '#include "povray_header.inc"'

    do i=lbound(particles,dim=1), ubound(particles,dim=1)
       do j=lbound(particles,dim=2), ubound(particles,dim=2)
          if (particles(i,j)%ID<=0) cycle
          call sphere(particles(i,j))
          nWritten = nWritten+1
       end do
    end do

    if (present(particles1)) then
       do i=lbound(particles1,dim=1), ubound(particles1,dim=1)
          do j=lbound(particles1,dim=2), ubound(particles1,dim=2)
             if (particles1(i,j)%ID<=0) cycle
             call sphere(particles1(i,j))
             nWritten = nWritten+1
          end do
       end do
    end if

    close(77)


    open(77,file=trim(filename)//'.dat')
    write(77,*) nWritten

    do i=lbound(particles,dim=1), ubound(particles,dim=1)
       do j=lbound(particles,dim=2), ubound(particles,dim=2)
          if (particles(i,j)%ID<=0) cycle
          call sphereDat(particles(i,j))
       end do
    end do

    if (present(particles1)) then
       do i=lbound(particles1,dim=1), ubound(particles1,dim=1)
          do j=lbound(particles1,dim=2), ubound(particles1,dim=2)
             if (particles1(i,j)%ID<=0) cycle
             call sphereDat(particles1(i,j))
          end do
       end do
    end if

    close(77)

  contains


    subroutine writeHeader
      open(78,file="povray_header.inc")
      ! write(78,'(A)') "#version 3.6;"
      write(78,'(A)') 'global_settings {  assumed_gamma 1.0 }'
      write(78,'(A)') '#default{ finish{ ambient 0.1 diffuse 0.9 }}'
      write(78,'(A)') '#include "colors.inc"'
      write(78,'(A)') '#include "textures.inc"'
      write(78,'(A)') '// camera ---------------------------------'
      write(78,'(A)') 'camera {location <33.0, 0., -33.0>'
      write(78,'(A)') '        look_at  < 0.0, 0.,   0.0>}'
      write(78,'(A)') '// sun ------------------------------------'
      write(78,'(A)') 'light_source{<300,100,-200> color White}'
      write(78,'(A)') '// macros ------------------------------------'
      write(78,'(A)') '#macro proton (rx,ry,rz)'
      write(78,'(A)') 'sphere{<rx,ry,rz>,1.00 texture{ pigment{color rgb<0,0,1>} finish{phong 1} } }'
      write(78,'(A)') '#end'
      write(78,'(A)') '#macro neutron (rx,ry,rz)'
      write(78,'(A)') 'sphere{<rx,ry,rz>,1.00 texture{ pigment{color rgb<1,0,0>} finish{phong 1} } }'
      write(78,'(A)') '#end'
      write(78,'(A)') '#macro antiNucleon (rx,ry,rz)'
      write(78,'(A)') 'sphere{<rx,ry,rz>,1.00 texture{ pigment{color rgb<1,1,1>} finish{phong 1} } }'
      write(78,'(A)') '#end'
      write(78,'(A)') '#macro pion (rx,ry,rz)'
      write(78,'(A)') 'sphere{<rx,ry,rz>,0.75 texture{ pigment{color rgb<0,1,0>} finish{phong 1} } }'
      write(78,'(A)') '#end'
      write(78,'(A)') '#macro Mes (rx,ry,rz)'
      write(78,'(A)') 'sphere{<rx,ry,rz>,0.75 texture{ pigment{color rgb<1,0,1>} finish{phong 1} } }'
      write(78,'(A)') '#end'
      write(78,'(A)') '#macro ExoMes (rx,ry,rz)'
      write(78,'(A)') 'sphere{<rx,ry,rz>,0.75 texture{ pigment{Candy_Cane} } }'
      write(78,'(A)') '#end'
      write(78,'(A)') '#macro Bar (rx,ry,rz)'
      write(78,'(A)') 'sphere{<rx,ry,rz>,1.25 texture{ pigment{color rgb<1,1,0>} finish{phong 1} } }'
      write(78,'(A)') '#end'
      write(78,'(A)') '#macro ExoBar (rx,ry,rz)'
      write(78,'(A)') 'sphere{<rx,ry,rz>,1.25 texture{ pigment{X_Gradient} } }'
      write(78,'(A)') '#end'
      close(78)
    end subroutine writeHeader


    subroutine writeScript
      use inputGeneral, only: numTimeSteps
      open(79,file="povray.sh")
      write(79,'(A,i4,A)') 'for ((  i = 1 ;  i <=',numTimeSteps,';  i++  ))'
      write(79,'(A)') 'do'
      write(79,'(A)') '  i_str=`printf %03i $i`'
      write(79,'(A)') '  povray Display=false Movie_Real_$i_str.pov'
      write(79,'(A)') '  povray Display=false Movie_Pert_$i_str.pov'
      write(79,'(A)') '  povray Display=false Movie_All_$i_str.pov'
      write(79,'(A)') 'done'
      write(79,'(A)') 'convert -delay 2 Movie_Real_* Movie_Real.gif'
      write(79,'(A)') 'convert -delay 2 Movie_Pert_* Movie_Pert.gif'
      write(79,'(A)') 'convert -delay 2 Movie_All_* Movie_All.gif'
      write(79,'(A)') 'rm Movie_*.png'
      close(79)
    end subroutine writeScript


    !**************************************************************************
    !****s* makePicture/sphere
    ! NAME
    ! subroutine sphere(iPart)
    ! INPUTS
    ! * type(particle) :: iPart
    ! PURPOSE
    ! For each particle a sphere is generated around its origin.
    !
    ! Color scheme:
    ! * proton : blue
    ! * neutron: red
    ! * pion  : green
    ! * other baryon : yellow (rgb<1,1,0>)
    ! * other meson : magenta (rgb<1,0,1>)
    ! * antiNucleon : black
    !
    ! Radius scheme:
    ! * nucleon : 1 fm
    ! * other baryon: 1.25 fm
    ! * meson: 0.75 fm
    !**************************************************************************
    subroutine sphere(iPart)
      use particleProperties, only: isCharmed, isStrange
      use IdTable, only: nucleon, Delta, pion

      type(particle),intent(in) :: iPart

      character(*), parameter :: T = '(A,"(",F8.3,",",F8.3,",",F8.3,")")'

      select case (iPart%ID)
      case (nucleon)
        if (iPart%antiparticle) then
          write(77,T) 'antiNucleon', iPart%position(1:3)
        else if (iPart%charge==1) then
          write(77,T) 'proton', iPart%position(1:3)
        else if (iPart%charge==0) then
          write(77,T) 'neutron', iPart%position(1:3)
        end if
      case (pion)
        write(77,T) 'pion', iPart%position(1:3)
      case (delta:pion-1)
        ! Baryons
        if (isCharmed(iPart%ID) .or. isStrange(iPart%ID)) then
          write(77,T) 'ExoBar', iPart%position(1:3)
        else
          write(77,T) 'Bar', iPart%position(1:3)
        end if
      case default
        if (isCharmed(iPart%ID) .or. isStrange(iPart%ID)) then
          write(77,T) 'ExoMes', iPart%position(1:3)
        else
          write(77,T) 'Mes', iPart%position(1:3)
        end if
      end select

    end subroutine sphere


    !**************************************************************************
    !****s* makePicture/sphereDat
    ! NAME
    ! subroutine sphereDat(iPart)
    ! INPUTS
    ! * type(particle) :: iPart
    ! PURPOSE
    ! For each particle the position and some additional info is written.
    !**************************************************************************
    subroutine sphereDat(iPart)
      type(particle), intent(in) :: iPart

      integer :: iType
      character(*), parameter :: T1 = '(1P,3e12.3,0P,i5,i3,f7.3)'

      iType = iPart%ID
      if (iPart%antiparticle) iType = -iPart%ID

      write(77,T1) iPart%position(1:3), iType, iPart%charge, iPart%scaleCS

    end subroutine sphereDat


  end subroutine makePicture


end module povray
