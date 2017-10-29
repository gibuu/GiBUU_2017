reset

# use svg terminal with mouse support (requires gnuplot 4.6)
set terminal svg enhanced mouse jsdir "http://gnuplot.sourceforge.net/demo_svg/" font "Helvetica,14"
set output 'Plot23.svg'

#set terminal postscript eps enhanced color 'Helvetica' 24
#set output 'Plot23.eps'
#set size 1,1
#set tics scale 2

set log xy
set yrange [1:300]
set xlabel "sqrt(s) [GeV]"
set ylabel "σ [mb]"

f = "XS.dat"

# Plotting as function of sqrt(s):

plot f u 2:($9<1?$5:1/0) w l lw 1 lt 3 not, \
     f u 2:($9<1?$6:1/0) w l lw 1 lt 3 not, \
     f u 2:($9>0?$7:1/0) w l lw 1 lt 4 not, \
     f u 2:($9>0?$8:1/0) w l lw 1 lt 4 not, \
     f u 2:3             w l lw 3 lt 1 t 'total', \
     f u 2:4             w l lw 3 lt 2 t 'elastic'
