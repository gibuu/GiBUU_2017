set terminal postscript eps enhanced color 'Helvetica' 14 dl 2

w = 3

mN=0.938
srts(x) = (2*mN*((mN**2+x**2)**0.5+mN))**0.5

xlab = 3.0
xmax = 4.0
set xrange [2:xmax]
set x2range [2:xmax]
set mxtics 5

set mytics 5

unset bars

set output 'XS_NN_1pi.eps'

set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

unset key

set multiplot layout 2,2

###################

set label 'pp {/Symbol \256} pn{/Symbol p^+}' at xlab,17

set format x ""
set x2tics

set ylabel "{/Symbol s}_{pp} [mb]"
set yrange [-5:21]

f = "pp_pnpi+.dat"
d = "LB_pp_pnpi+.dat"
d2 = "pp_pnpi+_subtracted.dat"

#f(x) = A * ((x-2.015)*5)**B * exp(-(C*((x-2.015)*5)**D+E*((x-2.015)*5)))

#fit f(x) d2 using 2:3 via A,B,C,D,E

plot \
     f u 1:($2+$8) w l lw w lc 1 t "total", \
     f u 1:3 w l lw w lc 2 t "Delta", \
     f u 1:4 w l lw w lc 3 t "I=1/2", \
     f u 1:5 w l lw w lc 4 t "I=3/2", \
     f u 1:8 w l lw w lc 7 t "BG Weil", \
     d u (srts($1)):2:3 w yerrorbars pt 7 lt 1 lc 0 lw w t 'data'

#     f u 1:6 w l lw w lc 5 t "BG Teis", \
#     f u 1:7 w l lw w lc 6 t "BG Buss", \
#     d2 u (srts($1)):3:4 w yerrorbars pt 6 lt 1 lc 1 t 'data subtr.'

###################

unset label
set label 'pp {/Symbol \256} pp{/Symbol p^0}' at xlab,4

unset ylabel
set format y ""
set yrange [-1:5]
set y2range [-1:5]
set y2tics

f = "pp_pppi0.dat"
d = "LB_pp_pppi0.dat"
d2 = "pp_pppi0_subtracted.dat"

replot

###################

unset label
set label 'pn {/Symbol \256} pn{/Symbol p^0}' at xlab,7

set format x "%g"
set xlabel "sqrt(s) [GeV]"
unset x2tics

set ylabel "{/Symbol s}_{pn} [mb]"
set format y "%g"
set yrange [0:8.8]
unset y2tics

f = "pn_pnpi0.dat"
d = "LB_pn_pnpi0.dat"
d2 = "pn_pnpi0_subtracted.dat"

replot

###################

unset label
set label 'pn {/Symbol \256} pp{/Symbol p^-}' at xlab,2.5

set key top right

unset ylabel
set format y ""
set yrange [-0.5:3.3]
set y2range [-0.5:3.3]
set y2tics

f = "pn_pppi-.dat"
d = "LB_pn_pppi-.dat"
d2 = "pn_pppi-_subtracted.dat"

replot

###################

unset multiplot


!./fixbb XS_NN_1pi.eps
