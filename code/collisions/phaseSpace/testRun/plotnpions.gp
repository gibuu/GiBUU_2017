reset

set xrange [0:2]

set xlabel "E_{lab} [GeV]"
set ylabel "prob"

#set log y
#set yrange [1e-3:]

plot "fort.141" u 1:($2/2) w l lw 6 t "2 pions"
replot "fort.142" u 1:($2/3) w l lw 6 t "3 pions"

set terminal postscript eps enhanced color 'Helvetica' 48
set output 'PlotNpions.eps'
set size 2,2
replot
set output
set terminal x11
