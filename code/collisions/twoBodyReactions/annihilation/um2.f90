SUBROUTINE USUN

               use inputGeneral, only: path_To_Input
               use output, only: Write_ReadingInput

!          ----- TO FILL THE ARRAY FOR ISOSCALAR FACTOR**2 --------
      LOGICAL TWO_NUCLEON,FITTING
      logical verbose
      COMMON /SWITCH/TWO_NUCLEON,FITTING
      COMMON /UN/ I2(2,4,2), IY2(2,4,2), U2(4,4,2),&
     &            I3(3,6,2), IY3(3,6,2), U3(4,6,2),&
     &            I4(4,9,2), IY4(4,9,2), U4(4,9,2),&
     &            I5(5,12,2),IY5(5,12,2),U5(4,12,2),&
     &            I6(6,16,2),IY6(6,16,2),U6(4,16,2),&
     & USUM(4,2,6), ISOS(2,6), IYSTR(2,6), NUMLIN(6)
      COMMON /VERB/ verbose
      SAVE /SWITCH/, /UN/, /VERB/

              if (TWO_NUCLEON) then
                open(10, file=trim(path_to_input)//'/annihilation/umd_tn.dat', STATUS='OLD')
              else
                open(10, file=trim(path_to_input)//'/annihilation/umd2.dat',   STATUS='OLD')
              end if

              if (verbose) call Write_ReadingInput('isoscalar factor tables',0)

              if (verbose) write(13,100)
  100  FORMAT(//6X,'SQUARES OF ISOSCALAR FACTORS TABLES: 2<N<6'/)

              do 20 N=2,6
              do 200 K=1,2
              do 2 I=1,4
           USUM(I,K,N)=0.
  2           CONTINUE
  200         CONTINUE
  20          CONTINUE

                           N=2
                    if (TWO_NUCLEON) then
                          NUMLIN(N)=2
                    else
                          NUMLIN(N)=4
                    end if
       do 12 K=1,2
         CALL HEAD(K,N)
         do 22 I=1,NUMLIN(N)
       READ (10,32)(I2(M,I,K),IY2(M,I,K),M=1,2),(U2(L,I,K),L=1,4)
  32   FORMAT(4I2,17X,4F11.6)
            do 42 I0=1,4
            USUM(I0,K,N)=USUM(I0,K,N)+U2(I0,I,K)
  42        CONTINUE
  22        CONTINUE
         CALL FOOT(K,N)
         do 72 I=1,NUMLIN(N)
            do 82 I0=1,4
            if (USUM(I0,K,N).LE.1.E-10) GOTO 82
            U2(I0,I,K)=U2(I0,I,K)/USUM(I0,K,N)
  82        CONTINUE
  72     CONTINUE
  12   CONTINUE



                           N=3
                    if (TWO_NUCLEON) then
                          NUMLIN(N)=4
                    else
                          NUMLIN(N)=6
                    end if
       do 13 K=1,2
         CALL HEAD(K,N)
         do 23 I=1,NUMLIN(N)
       READ (10,33)(I3(M,I,K),IY3(M,I,K),M=1,3),(U3(L,I,K),L=1,4)
  33   FORMAT(6I2,13X,4F11.6)
            do 43 I0=1,4
            USUM(I0,K,N)=USUM(I0,K,N)+U3(I0,I,K)
  43        CONTINUE
  23        CONTINUE
         CALL FOOT(K,N)
         do 73 I=1,NUMLIN(N)
            do 83 I0=1,4
            if (USUM(I0,K,N).LE.1.E-10) GOTO 83
            U3(I0,I,K)=U3(I0,I,K)/USUM(I0,K,N)
  83        CONTINUE
  73     CONTINUE
  13   CONTINUE



                           N=4
                    if (TWO_NUCLEON) then
                          NUMLIN(N)=6
                    else
                          NUMLIN(N)=9
                    end if
       do 14 K=1,2
         CALL HEAD(K,N)
         do 24 I=1,NUMLIN(N)
       READ (10,34)(I4(M,I,K),IY4(M,I,K),M=1,4),(U4(L,I,K),L=1,4)
  34   FORMAT(8I2,9X,4F11.6)
            do 44 I0=1,4
            USUM(I0,K,N)=USUM(I0,K,N)+U4(I0,I,K)
  44        CONTINUE
  24        CONTINUE
         CALL FOOT(K,N)
         do 74 I=1,NUMLIN(N)
            do 84 I0=1,4
            if (USUM(I0,K,N).LE.1.E-10) GOTO 84
            U4(I0,I,K)=U4(I0,I,K)/USUM(I0,K,N)
  84        CONTINUE
  74     CONTINUE
  14   CONTINUE



                           N=5
                    if (TWO_NUCLEON) then
                          NUMLIN(N)=9
                    else
                          NUMLIN(N)=12
                    end if
       do 15 K=1,2
         CALL HEAD(K,N)
         do 25 I=1,NUMLIN(N)
       READ (10,35)(I5(M,I,K),IY5(M,I,K),M=1,5),(U5(L,I,K),L=1,4)
  35   FORMAT(10I2,5X,4F11.6)
            do 45 I0=1,4
            USUM(I0,K,N)=USUM(I0,K,N)+U5(I0,I,K)
  45        CONTINUE
  25        CONTINUE
         CALL FOOT(K,N)
         do 75 I=1,NUMLIN(N)
            do 85 I0=1,4
            if (USUM(I0,K,N).LE.1.E-10) GOTO 85
            U5(I0,I,K)=U5(I0,I,K)/USUM(I0,K,N)
  85        CONTINUE
  75     CONTINUE
  15   CONTINUE



                           N=6
                    if (TWO_NUCLEON) then
                          NUMLIN(N)=12
                    else
                          NUMLIN(N)=16
                    end if
       do 16 K=1,2
         CALL HEAD(K,N)
         do 26 I=1,NUMLIN(N)
       READ (10,36)(I6(M,I,K),IY6(M,I,K),M=1,6),(U6(L,I,K),L=1,4)
  36   FORMAT(12I2,1X,4F11.6)
            do 46 I0=1,4
            USUM(I0,K,N)=USUM(I0,K,N)+U6(I0,I,K)
  46        CONTINUE
  26        CONTINUE
         CALL FOOT(K,N)
         do 76 I=1,NUMLIN(N)
            do 86 I0=1,4
            if (USUM(I0,K,N).LE.1.E-10) GOTO 86
            U6(I0,I,K)=U6(I0,I,K)/USUM(I0,K,N)
  86        CONTINUE
  76     CONTINUE
  16   CONTINUE

         if (verbose) call Write_ReadingInput('isoscalar factor tables',1)
         close(10)

         RETURN
         END



            SUBROUTINE HEAD(K,N)
            logical verbose
      COMMON /UN/ I2(2,4,2), IY2(2,4,2), U2(4,4,2),&
     &            I3(3,6,2), IY3(3,6,2), U3(4,6,2),&
     &            I4(4,9,2), IY4(4,9,2), U4(4,9,2),&
     &            I5(5,12,2),IY5(5,12,2),U5(4,12,2),&
     &            I6(6,16,2),IY6(6,16,2),U6(4,16,2),&
     & USUM(4,2,6), ISOS(2,6), IYSTR(2,6), NUMLIN(6)
      COMMON /VERB/ verbose
      SAVE /UN/, /VERB/
            READ (10,100)
 100        FORMAT(A70)
            if (verbose) write(13,100)
            READ (10,101)ISOS(K,N),IYSTR(K,N)
 101        FORMAT(15x,I2,5x,I2)
            if (verbose) write(13,101)ISOS(K,N),IYSTR(K,N)
            READ (10,100)
        RETURN
        END
                SUBROUTINE FOOT(K,N)
                logical verbose
      COMMON /UN/ I2(2,4,2), IY2(2,4,2), U2(4,4,2),&
     &            I3(3,6,2), IY3(3,6,2), U3(4,6,2),&
     &            I4(4,9,2), IY4(4,9,2), U4(4,9,2),&
     &            I5(5,12,2),IY5(5,12,2),U5(4,12,2),&
     &            I6(6,16,2),IY6(6,16,2),U6(4,16,2),&
     & USUM(4,2,6), ISOS(2,6), IYSTR(2,6), NUMLIN(6)
      COMMON /VERB/ verbose
      SAVE /UN/, /VERB/
            if (verbose) write(13,190)(USUM(I,K,N),I=1,4)
 190        FORMAT(27x,4F8.2)
              do 200 L=1,3
               READ (10,201)
201           FORMAT(A10,61x)
 200          CONTINUE
            RETURN
            END
