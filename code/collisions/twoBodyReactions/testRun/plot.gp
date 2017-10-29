set macro

#set terminal postscript eps enhanced color 'Helvetica' 24
#set output 'Plot.eps'

set terminal svg enhanced mouse size 900,800 jsdir "http://gnuplot.sourceforge.net/demo_svg/" font "Helvetica,14"

set log xy
set mxtics 10
set mytics 10
set ylabel "σ [mb]"
set yrange [1:400]

unset bars

f = 'XS.dat'

#################################################
# data files

DT = "dataTotal.dat"
DE = "dataElast.dat"

# extract projectile and target masses from data files (if present)
mp = `! if [ -s @DT ] ; then cat @DT | head -10 | tail -1 | awk '{print $1}' ; else echo "0" ; fi`
mt = `! if [ -s @DT ] ; then cat @DT | head -10 | tail -1 | awk '{print $2}' ; else echo "0" ; fi`

srts(p)= (((p**2+mp**2)**0.5 + mt)**2 - p**2)**0.5

#################################################

set output 'plot.svg'
set xlabel "sqrt(s) [GeV]"

plot \
     DT u (srts($2)):5:(srts($3)):(srts($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 7 ps 0.5 t 'data (total)', \
     DE u (srts($2)):5:(srts($3)):(srts($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 6 ps 0.5 t 'data (elastic)', \
     f u 2:($9<1?$5:1/0) w l lw 1 lt 3 t '', \
     f u 2:($9<1?$6:1/0) w l lw 1 lt 3 t '', \
     f u 2:($9>0?$7:1/0) w l lw 1 lt 4 t '', \
     f u 2:($9>0?$8:1/0) w l lw 1 lt 4 t '', \
     f u 2:3 w l lw 3 lt 1 t 'total', \
     f u 2:4 w l lw 3 lt 2 t 'elastic'

#################################################

set output 'plot_plab.svg'
set xlabel "p_{lab} [GeV]"

plot \
     DT u 2:5:3:4:($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 7 ps 0.5 t 'data (total)', \
     DE u 2:5:3:4:($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 6 ps 0.5 t 'data (elastic)', \
     f u 1:($9<1?$5:1/0) w l lw 1 lt 3 t '', \
     f u 1:($9<1?$6:1/0) w l lw 1 lt 3 t '', \
     f u 1:($9>0?$7:1/0) w l lw 1 lt 4 t '', \
     f u 1:($9>0?$8:1/0) w l lw 1 lt 4 t '', \
     f u 1:3 w l lw 3 lt 1 t 'total', \
     f u 1:4 w l lw 3 lt 2 t 'elastic'

#################################################
