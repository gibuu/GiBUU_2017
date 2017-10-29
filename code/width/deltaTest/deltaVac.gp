set terminal postscript eps enhanced color size 2.5,2.5 dl 2

set xlabel 'm_{/Symbol D} [GeV]'
set xrange [1:2.2]
set mxtics 2

set style data lines

f1 = 'DeltaVacWidth.1.dat'
f2 = 'DeltaVacWidth.2.dat'
f3 = 'DeltaVacWidth.3.dat'
f4 = 'DeltaVacWidth.4.dat'
f5 = 'DeltaVacWidth.5.dat'

set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0

w = 3

#######################################

set output 'deltaVacSF.eps'

set ylabel 'A_{/Symbol D} [GeV^{-2}]'
set mytics 2

set key box

set arrow 1 from 1.232,0 to 1.232,0.14 nohead lt 0

plot f1 u 1:3 lw w t 'Manley', \
     f2 u 1:3 lw w t 'Dmitriev' , \
     f3 u 1:3 lw w t 'Moniz' , \
     f4 u 1:3 lw w t 'Verwest', \
     f5 u 1:3 lw w t 'Bass'
     
!fixbb deltaVacSF.eps
!epstopdf deltaVacSF.eps

#######################################

set output 'deltaVacWidth.eps'

set ylabel '{/Symbol G_D} [GeV]'

set yrange [1E-3:1E1]

set log y
set format y '10^{%L}'
set mytics 10

set key top left

set arrow 1 from 1.232,1E-3 to 1.232,0.12 nohead lt 0
set arrow 2 from 1,0.12 to 1.232,0.12 nohead lt 0

plot f1 u 1:2 lw w t 'Manley', \
     f2 u 1:2 lw w t 'Dmitriev' , \
     f3 u 1:2 lw w t 'Moniz' , \
     f4 u 1:2 lw w t 'Verwest', \
     f5 u 1:2 lw w t 'Bass'

!fixbb deltaVacWidth.eps
!epstopdf deltaVacWidth.eps

#######################################

!rm *.eps
