!******************************************************************************
!****m* /ThreeBarMes
! NAME
! module ThreeBarMes
!
! PURPOSE
! This module implements the back reaction 3 -> 2 of baryon meson collisions
!
! INPUTS
! (none)
!******************************************************************************
module ThreeBarMes

  use callstack, only: traceBack

  implicit none

  private

  public :: doThreeBarMes

contains

  subroutine doThreeBarMes(time)

    real, intent(in) :: time

    call doNpipi(time)

  end subroutine doThreeBarMes

  subroutine doNpipi(time)

    use idTable, only: nucleon, pion, phi, omegaMeson
    use DecayChannels, only: tDecay3body
    use particleDefinition
    use VolumeElements, only: VolumeElements_boxSize, VolumeElements_3Body
    use particleProperties, only: hadron

    real, intent(in) :: time

    integer :: qNuc, qPi1, qPi2, iAnti
    type(tDecay3Body) :: decay
    integer :: isSame

    logical :: doInit
    integer, dimension(1:3) :: iEns, iInd
    type(particle), POINTER :: Part1, Part2, Part3
    real :: srts0, srts, scaleFak
    real, dimension(1:3) :: betaToCM

    srts0 = min(1.95, &
         hadron(omegaMeson)%mass + hadron(nucleon)%mass &
         &  - 3*hadron(omegaMeson)%width, &
         hadron(phi)%mass + hadron(nucleon)%mass &
         &  - 3*hadron(phi)%width )

    decay%ID = (/ nucleon, pion, pion /)

    do qNuc=0,1
       do iAnti=0,1
          select case (iAnti)
          case (0)
             decay%isAnti(1) = .false.
             decay%charge(1) = qNuc
          case (1)
             decay%isAnti(1) = .true.
             decay%charge(1) = -qNuc
          end select

          do qPi1 = -1,1
             decay%charge(2) = qPi1
             do qPi2 = qPi1,1
                decay%charge(3) = qPi2

                if (abs(sum(decay%charge)) == 3) cycle

                if (qPi1 /= qPi2) then
                   isSame = 0
                else
                   isSame = 6 ! pions are equal
                end if

                call doMain

             end do
          end do
       end do
    end do

  contains

    subroutine doMain

      use constants, only: pi, hbarc, mPi
      use mediumDefinition
      use preEventDefinition, only: preEvent
      use finalStateModule, only: massAss
      use lorentzTrafo, only: lorentzCalcBeta !, lorentz
!      use output, only: WriteParticle
      use threeMeson, only: calcPhi3
      use inputGeneral, only: delta_T, numEnsembles
      use random, only: rn
      use mesonWidth, only: calcSpectralIntegral!, getSpectralIntegral
      use twoBodyTools, only: pCM

      use pionNucleon, only: pionNuc
      use phiNucleon, only: phiNuc
      use omegaNucleon, only: omegaNuc


      integer :: iq1,q1,q2
      type(particle), dimension(1:4) :: partIn
      type(particle), dimension(1:2) :: partOut
      type(medium) :: mediumATcoll
      real, dimension(0:3) :: momLRF,momCalc
      real, dimension(1:3) :: betaToLRF
      logical :: flagOk
      real :: sigmaTot, sigmaElast
      real, dimension(1:20) :: sigmaArr
      type(preEvent),dimension(1:4) :: preOut
      real :: preFak, Prob31, XS

      ! cf. master_2body.f90;  LRF: defined by density
      ! momLRF:    sum of all momenta in LRF
      ! momCalc:   sum of all momenta in calc frame
      ! betaToLRF: beta of LRF in calc frame
      ! betaToCM:  beta of CM in calc frame


      call SetToDefault(partIn) ! just dummy values!

      doInit = .true.
      tripleLoop: do
         if (.not.VolumeElements_3Body(decay, isSame, &
              doInit, iEns,iInd, Part1,Part2,Part3, scaleFak, 2)) exit

         srts = sqrtS(Part1,Part2,Part3)
         if (srts < srts0) cycle tripleLoop

         momCalc = Part1%momentum+Part2%momentum+Part3%momentum
         betaToCM = lorentzCalcBeta(momCalc)

         ! we assume: LRF = calc frame
         momLRF = momCalc
         betaToLRF = 0.

         preFak = delta_T/(numEnsembles*VolumeElements_boxSize())**2 &
              * 0.1*hbarc**3 * 6.0 * scaleFak &
              / ( 2 * Part1%momentum(0)*Part2%momentum(0)*Part3%momentum(0) &
              * 2.0 * calcPhi3(srts, mPi) * 4*pi)

         do iq1 = 0,1
            q1 = iq1
            if (iAnti==1) q1 = -iq1
            q2 = sum(decay%charge) - q1

            partOut(1)%ID = nucleon
            partOut(1)%charge = q1
            partOut(1)%antiparticle = (iAnti==1)

            partOut(2)%charge = q2
            partOut(2)%antiparticle = .false.

            select case (q2)
            case (0) !===== pi0, omega, phi =====

!               write(*,*) '3BM: ',iAnti,qNuc,qPi1,qPi2, srts, scaleFak, &
!                    calcSpectralIntegral(omegaMeson,srts-hadron(nucleon)%mass),&
!                    calcSpectralIntegral(phi,srts-hadron(nucleon)%mass)

               !--- pi0

               if (srts>1.95) then
                  partOut(2)%ID = pion
                  call massAss(srts, mediumAtColl, partIn(1:2), partOut(1:2), &
                       betaToLRF,betaToCM,0,flagOk)

                  if (flagOk) then

!                     call writeParticle(6,99,partOut)

                     ! we use partOut as the incoming (!) particles
                     call pionNuc(srts, partOut, mediumAtColl, momLRF, &
                          preOut, sigmaTot, sigmaElast, .false., sigmaArr)

                     XS = sigmaArr(19)
                     Prob31 = preFak * XS &
                          * pCM(srts, partOut(1)%mass, partOut(2)%mass)

                     if (Prob31 > 1.0) then
                        write(*,*) 'P > 1: pi pi N -> pi N  :',Prob31
                     end if
                     if (Prob31 > rn()) then
                        call doColl(Part1,Part2,Part3, partOut(1:2), time)
                        cycle tripleLoop
                     end if


                  end if
               end if

               !--- omega

               ! we cut it arbitrarily:
               if (srts > hadron(omegaMeson)%mass + hadron(nucleon)%mass &
                    - 3*hadron(omegaMeson)%width) then

                  partOut(2)%ID = omegaMeson

                  call massAss(srts, mediumAtColl, partIn(1:2), partOut(1:2), &
                       betaToLRF,betaToCM,0,flagOk)

                  if (flagOk) then

!                     call writeParticle(6,99,partOut)

                     ! we use partOut as the incoming (!) particles
                     call omegaNuc(srts, partOut, mediumAtColl, momLRF, &
                          preOut(1:3), sigmaTot, sigmaElast, &
                          .false., 0.0, & ! <-- TO BE DONE
                          .false., &
                          1.0, &  ! <-- TO BE DONE
                          sigmaArr(1:6))

                     XS = sigmaArr(3)
                     Prob31 = preFak * XS &
                          * pCM(srts, partOut(1)%mass, partOut(2)%mass) &
                          * calcSpectralIntegral(omegaMeson,srts-hadron(nucleon)%mass)

                     if (Prob31 > 1.0) then
                        write(*,*) 'P > 1: pi pi N -> omega N  :',Prob31
                     end if
                     if (Prob31 > rn()) then
                        call doColl(Part1,Part2,Part3, partOut(1:2), time)
                        cycle tripleLoop
                     end if

                  end if
               end if

               !--- phi

               ! we cut it arbitrarily:
               if (srts > hadron(phi)%mass + hadron(nucleon)%mass &
                    - 3*hadron(phi)%width) then

                  partOut(2)%ID = phi

                  call massAss(srts, mediumAtColl, partIn(1:2), partOut(1:2), &
                       betaToLRF,betaToCM,0,flagOk)

                  if (flagOk) then

!                     call writeParticle(6,99,partOut)

                     ! we use partOut as the incoming (!) particles
                     call phiNuc(srts, partOut, mediumAtColl, &
                          preOut(1:3), sigmaTot, sigmaElast, &
                          .false., 0.0, & ! <-- TO BE DONE
                          .false., sigmaArr(1:3))

                     XS = sigmaArr(3)
                     Prob31 = preFak * XS &
                          * pCM(srts, partOut(1)%mass, partOut(2)%mass) &
                          * calcSpectralIntegral(phi,srts-hadron(nucleon)%mass)

                     if (Prob31 > 1.0) then
                        write(*,*) 'P > 1: pi pi N -> phi N  :',Prob31
                     end if
                     if (Prob31 > rn()) then
                        call doColl(Part1,Part2,Part3, partOut(1:2), time)
                        cycle tripleLoop
                     end if
                  end if
               end if

               !--- finished


!               stop

            case (-1,1) !===== pi-,pi+ =====

               if (srts>1.95) then
                  partOut(2)%ID = pion
                  call massAss(srts, mediumAtColl, partIn(1:2), partOut(1:2), &
                       betaToLRF,betaToCM,0,flagOk)

                  if (flagOk) then

!                     call writeParticle(6,99,partOut)

                     ! we use partOut as the incoming (!) particles
                     call pionNuc(srts, partOut, mediumAtColl, momLRF, &
                          preOut, sigmaTot, sigmaElast, .false., sigmaArr)


                     XS = sigmaArr(19)
                     Prob31 = preFak * XS &
                          * pCM(srts, partOut(1)%mass, partOut(2)%mass)

                     if (Prob31 > 1.0) then
                        write(*,*) 'P > 1: pi pi N -> pi N  :',Prob31
                     end if
                     if (Prob31 > rn()) then
                        call doColl(Part1,Part2,Part3, partOut(1:2), time)
                        cycle tripleLoop
                     end if


                  end if
               end if


            case default
               !... do nothing

            end select
         end do


      end do tripleLoop

    end subroutine doMain

    subroutine doColl(part1,part2,part3, partOut, time)

      use collisionNumbering, only: real_numbering, real_firstnumbering, &
           ReportEventNumber
      use history, only: setHistory
      use twoBodyStatistics, only: rate
      use particleDefinition
!      use output, only: WriteParticle
      use collisionTools, only: finalcheck
      use lorentzTrafo, only: lorentz

      type(particle), POINTER, intent(inout) :: part1, part2, part3
      type(particle), dimension(1:2), intent(inout) :: partOut
      real, intent(in) :: time

      integer :: number
      logical :: flag


      !===== Set new particles =====

      partOut(1)%position = (part1%position+part2%position+part3%position)/3
      partOut(2)%position = partOut(1)%position

      partOut(1:2)%number = 0 ! first reset it, then set it again
      call setNumber(partOut)

      call setHistory(part1,part2,part3, partOut)

      ! no formation time effects:
      partOut%lastCollisionTime = time
      partOut%productionTime    = time
      partOut%formationTime     = time
      !    partOut%scaleCS=1.
      !    partOut%in_Formation=.false.

      ! particles are not perturbative:
      !    partOut%perturbative = .false.
      !    partOut%perWeight = 1.0

      ! particles are 'on-shell':
      !    partOut%offshellParameter = 0.

      ! now the tricky part starts:
      number=real_numbering()
      partOut%event(1)=number
      partOut%event(2)=number

      partOut%firstEvent = real_firstnumbering()

      call lorentz(-betaToCM,partOut(1)%momentum)
      call lorentz(-betaToCM,partOut(2)%momentum)

      !===== Set some global counters =====

      call reportEventNumber( (/Part1,Part2,Part3/), PartOut, PartOut(1)%event, time, 3111)

      call rate( (/Part1,Part2,Part3/), PartOut, time)

      flag = finalCheck((/Part1,Part2,Part3/), PartOut)
      if (.not.flag) call Traceback('finalCheck')

      !===== Set particle 1 & 2, Delete particle 3 =====

!      write(*,*) '3BM: ',iAnti,qNuc,qPi1,qPi2, srts, scaleFak
!      call writeParticle(6,99,partOut)

      part1 = partOut(1)
      part2 = partOut(2)
      call setToDefault(part3)


    end subroutine doColl

  end subroutine doNpipi

end module ThreeBarMes
