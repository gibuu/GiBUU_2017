
set terminal postscript eps enhanced dashed

set output 'DeltaDalitzFF.eps'

f = 'DeltaDalitzFF.dat'

set xlabel 'm_{ee} [GeV]'

set xrange [0:1.5]
set mxtics 2

set ylabel '|G_M|^2'

set log y
set yrange [4e-1:4e2]
set format y "10^{%L}"
set mytics 10

set key top right box opaque

w = 3

set style data lines

plot \
     3.029**2 lt 0, \
     f u 1:4  lt 1 lc 1 lw w t 'W/I 1.23', \
     f u 1:5  lt 2 lc 1 lw w t 'W/I 1.43', \
     f u 1:6  lt 3 lc 1 lw w t 'W/I 1.63', \
     f u 1:7  lt 3 lc 1 lw w t 'W/I 1.83', \
     f u 1:8  lt 4 lc 1 lw w t 'W/I 2.03', \
     f u 1:9  lt 1 lc 2 lw w t 'R/P 1.23', \
     f u 1:10 lt 2 lc 2 lw w t 'R/P 1.43', \
     f u 1:11 lt 3 lc 2 lw w t 'R/P 1.63', \
     f u 1:12 lt 3 lc 2 lw w t 'R/P 1.83', \
     f u 1:13 lt 4 lc 2 lw w t 'R/P 2.03' 

#     f u 1:2  lt 1 lw w t 'VMD', \
#     f u 1:3  lt 2 lw w t 'Dipole', \