
set terminal postscript eps enhanced dashed color

set output 'DeltaDalitzFF_Wdep.eps'

f = 'DeltaDalitzFF_Wdep.dat'

set xlabel 'W [GeV]'

set xrange [1:2.1]
set mxtics 2

set ylabel '|G_M|^2'

#set log y
set yrange [0:11]
#set format y "10^{%L}"
#set mytics 10

set key top right box opaque

w = 3

set style data lines

set arrow from 1.232,0 to 1.232,11 nohead lt 0

plot \
     f u 1:2 lw w t 'Iachello', \
     f u 1:3 lw w t 'Ramalho', \
     f u 1:($4*3**2) lw w t 'Compton'
