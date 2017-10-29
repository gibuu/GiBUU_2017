program main

  use inputGeneral, only: readInputGeneral
  use particleProperties, only: initParticleProperties
  use particleDefinition
  use history, only: setHistory, history_getParents, history_getGeneration
  use IdTable, only: nucleon, pion, omegaMeson

  implicit none

  integer, dimension(1:3) :: parents
  type(particle) :: p1, p2, p3, part(1)
!   integer :: hist

  call initParticleProperties
  call readInputGeneral


  p1%ID = nucleon
  p1%history = 1000000

  p2%ID = 0
  p2%history = 0

  p3%ID = 0
  p3%history = 0

  call setHistory (p1, p2, p3, part)
  write(*,*)  "history    =", part%history
  parents = history_getParents(part(1)%history)
  write(*,*)  "parents    =", parents 
  write(*,*)  "generation =", history_getGeneration(part(1)%history)
!   write(*,*)  "1Body =", history_1Body(hist)
!   write(*,*)  "2Body =", history_2Body(hist)
!   write(*,*)  "3Body =", history_3Body(hist)
  write(*,*)


  p1%ID = nucleon
  p1%history = 1000000

  p2%ID = omegaMeson
  p2%history = 2000121

  p3%ID = 0
  p3%history = 0

  call setHistory (p1, p2, p3, part)
  write(*,*)  "history    =", part%history
  parents = history_getParents(part(1)%history)
  write(*,*)  "parents    =", parents 
  write(*,*)  "generation =", history_getGeneration(part(1)%history)
!   write(*,*)  "1Body =", history_1Body(hist)
!   write(*,*)  "2Body =", history_2Body(hist)
!   write(*,*)  "3Body =", history_3Body(hist)
  write(*,*)


  p1%ID = nucleon
  p1%history = 1000000

  p2%ID = pion
  p2%history = 2000121

  p3%ID = nucleon
  p3%history = 13120

  call setHistory (p1, p2, p3, part)
  write(*,*)  "history    =", part%history
  parents = history_getParents(part(1)%history)
  write(*,*)  "parents    =", parents 
  write(*,*)  "generation =", history_getGeneration(part(1)%history)
!   write(*,*)  "1Body =", history_1Body(hist)
!   write(*,*)  "2Body =", history_2Body(hist)
!   write(*,*)  "3Body =", history_3Body(hist)
  write(*,*)


end program main