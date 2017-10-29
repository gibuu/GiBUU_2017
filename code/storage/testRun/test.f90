program main

  use PartInfoList_nLead


  real :: r
  integer :: i
  logical :: flag

  call PartInfoList_nLead_PUT(3,0.33)
  call PartInfoList_nLead_PUT(5,0.55)
  call PartInfoList_nLead_PUT(7,0.77)
  call PartInfoList_nLead_PUT(8,0.88)


  call PartInfoList_nLead_PUT(1,0.11)
  call PartInfoList_nLead_PUT(6,0.66)
  call PartInfoList_nLead_PUT(7,0.314)


  do i=1,10
     flag = PartInfoList_nLead_GET(i,r)
     write(*,*) i,flag, r
  enddo

  write(*,*) ' ========'
  call PartInfoList_nLead_ZERO
  call PartInfoList_nLead_PUT(7,0.77)
  
  do i=1,10
     flag = PartInfoList_nLead_GET(i,r)
     write(*,*) i,flag, r
  enddo


end program main
