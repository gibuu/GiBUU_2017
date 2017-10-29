
set terminal postscript eps enhanced dashed color

set xlabel 'm_{ee} [GeV]'

set mxtics 2

set ylabel 'd{/Symbol G}/dm_{ee}'

set log y
set yrange [1e-8:1e-2]
set format y "10^{%L}"
set mytics 10

w = 3

set style data lines

#######################################

set output 'DeltaDalitz.eps'

f = 'dGamma_dM_DeltaDalitz.dat'

plot \
     f u 1:7  lt 1 lc 1 lw w t 'Krivoruchenko (W = 1.23)', \
     f u 1:2  lt 2 lc 1 lw w t 'Wolf (W = 1.23)', \
     f u 1:3  lt 2 lc 2 lw w not, \
     f u 1:4  lt 2 lc 3 lw w not, \
     f u 1:5  lt 2 lc 4 lw w not, \
     f u 1:6  lt 2 lc 5 lw w not, \
     f u 1:12 lt 3 lc 1 lw w t 'Ernst (W = 1.23)', \
     f u 1:13 lt 3 lc 2 lw w not, \
     f u 1:14 lt 3 lc 3 lw w not, \
     f u 1:15 lt 3 lc 4 lw w not, \
     f u 1:16 lt 3 lc 5 lw w not, \
     f u 1:8  lt 1 lc 2 lw w t 'W = 1.43', \
     f u 1:9  lt 1 lc 3 lw w t 'W = 1.63', \
     f u 1:10 lt 1 lc 4 lw w t 'W = 1.83', \
     f u 1:11 lt 1 lc 5 lw w t 'W = 2.03'

#######################################

set output 'DeltaDalitz_rhoN.eps'

f = 'dGamma_dM_rhoN.dat'

plot \
    f u 1:3 lt 1 lw w t 'W = 1.23', \
    f u 1:4 lt 2 lw w t 'W = 1.43', \
    f u 1:5 lt 3 lw w t 'W = 1.63', \
    f u 1:6 lt 4 lw w t 'W = 1.83', \
    f u 1:7 lt 5 lw w t 'W = 2.03'

#######################################

!epstopdf DeltaDalitz.eps
!epstopdf DeltaDalitz_rhoN.eps

!rm *.eps
