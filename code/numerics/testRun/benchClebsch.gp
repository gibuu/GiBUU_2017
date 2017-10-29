
set terminal pdf enhanced

set output 'benchClebsch.pdf'

set style data lines

set key top left box

set log y
set format y "10^{%L}"

set mxtics 2
set mytics 10

set xlabel 'j_{max}'
set ylabel 't [s]'

w = 3

f = 'bench.dat'

plot \
    f u 1:2 lw w t 'GSL', \
    f u 1:3 lw w t 'GiBUU'
