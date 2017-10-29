set macro

set terminal postscript eps enhanced color dl 2 "Helvetica" 18

set style data lines
w = 3

set xlabel 'sqrt(s) [GeV]'
set ylabel '{/Symbol s} [mb]'

################################################################################

set output 'NNeta.eps' 

set title 'Res. Model: p p {/Symbol \256} p p {/Symbol h} X'

f1 = 'pp_eta.dat'
f2 = 'pp_pi_eta.dat'
d = 'LB_eta_excl.dat'

set xrange [2:4.5]
set yrange [1e-2:1e1]
set log y
set key top left box

mN=0.938
srts(x) = (2*mN*((mN**2+x**2)**0.5+mN))**0.5

plot f1 u 1:($2+$3) lt 1 lc 2 lw w t 'pp{/Symbol h}', \
     f2 u 1:2       lt 2 lc 2 lw w t 'NN{/Symbol ph}', \
     d u (srts($1)):2:3 w yerrorbars pt 7 lt 1 lc 2 lw w t 'data (LB,excl)'

#     f1 u 1:2       lw w t 'I=1/2', \
#     f1 u 1:3       lw w t 'I=3/2', \
#     f1 u 1:4       lw w t 'S11(1535)', \



################################################################################

set output 'NNrho.eps' 

set title 'p p {/Symbol \256} p R {/Symbol \256} p p {/Symbol r^0}'

f = 'pp_rho.dat'

set xrange [2:5.5]
set yrange [0:0.225]
unset log y
set mytics 2

plot f u 1:($2+$3) lw w lt 1 lc 1 t 'total', \
     f u 1:2       lw w lt 1 lc 2 t 'N^*', \
     f u 1:3       lw w lt 1 lc 3 t '{/Symbol D}^*', \
     "rusta_rho_excl.png.dat" u 1:2     w points     pt 6 lt 1 lc 0 lw w not, \
     "LB_rho_excl.dat" u (srts($1)):2:3 w yerrorbars pt 7 lt 1 lc 0 lw w t 'data'

#     f u 1:9       lw w lt 2 lc 2 t 'D13(1520)', \
#     f u 1:6       lw w lt 3 lc 2 t 'S11(1535)', \
#     f u 1:7       lw w lt 4 lc 2 t 'S11(1650)', \
#     f u 1:18      lw w lt 5 lc 2 t 'F15(1680)', \
#     f u 1:16      lw w lt 6 lc 2 t 'P13(1720)', \
#     f u 1:13      lw w lt 7 lc 2 t 'G17(2190)', \
#     f u 1:17      lw w lt 8 lc 2 t 'P13(1900)', \

#     f u 1:21      lw w lt 2 lc 3 t 'S31(1620)', \
#     f u 1:23      lw w lt 3 lc 3 t 'D33(1700)', \
#     f u 1:28      lw w lt 4 lc 3 t 'P31(1910)', \
#     f u 1:32      lw w lt 5 lc 3 t 'F35(1905)'

#     f u 1:($6+$7+$9+$16+$18+$21+$23+$28+$32) lw 1 lt 1 lc 4 t 'tot', \
#     f u 1:($6+$7+$9+$16+$18) lw 1 lt 1 lc 4 t 'I=1/2', \
#     f u 1:($21+$23+$28+$32) lw 1 lt 1 lc 4 t 'I=3/2', \

################################################################################

set output 'NNomega.eps' 

set title 'p p {/Symbol \256} p p {/Symbol w} X'

F1 = 'pp_omega.dat'
F2 = 'pp_pi_omega.dat'
FF = 'omega.dat'

!paste @F2 @F1 > @FF

set xrange [2.5:5]
set yrange [0:0.8]
unset log y
set mytics 2

set key top left

sqrts0 = 2.658
f(x) = 2.5*((x/sqrts0)**2-1.)**1.47*((x/sqrts0)**2)**-1.11

mN=0.938
srts(x) = (2*mN*((mN**2+x**2)**0.5+mN))**0.5

plot \
     FF u 1:($2+$17+$21) lw w lt 1 lc 1 t 'total', \
     FF u 1:($17+$21) lw w lt 1 lc 2 t 'p p {/Symbol w}', \
     FF u 1:21        lw w lt 2 lc 2 t 'N P_{13}(1900)', \
     FF u 1:17        lw w lt 3 lc 2 t 'N G_{17}(2190)', \
     FF u 1:2         lw w lt 1 lc 3 t 'p p {/Symbol p w}', \
     FF u 1:3         lw w lt 2 lc 3 t '{/Symbol D} P_{13}(1900)', \
     FF u 1:4         lw w lt 3 lc 3 t '{/Symbol D} G_{17}(2190)', \
     "ppomega.png.dat" u (sqrts0+$1/1E3):($2/1E3) w points lt 1 lc 2 pt 5 lw w not, \
     "rusta_omega.png.dat" u 1:2 w points pt 6 lw w lt 1 lc 2 not, \
     "LaBoer.dat" u (srts($1)):2:3 w yerrorbars pt 7 lt 1 lc 2 lw w t 'data ({/Symbol w}pp)', \
     "HADES.dat" u 1:2:3 w yerrorbars pt 6 lw w lt 1 lc 1 t 'HADES({/Symbol w}X)'

#      f(x) lw w lt 3 lc 1 t'Sibirtsev({/Symbol w}X)'

#     f u 1:2       lw w lt 1 lc 2 t 'I=1/2', \
#     f u 1:3       lw w lt 1 lc 3 t 'I=3/2', \

