set macro

set terminal postscript eps enhanced color dl 2 "Helvetica" 18

set yrange [0:60]
set mytics 2

set ylabel "{/Symbol s}_{pp} [mb]"

set style data lines
w = 3
w3 = 5
unset bars

d1 = 'PDG_pp_total.dat'
d2 = 'PDG_pp_elastic.dat'

f1 = 'pp_NR_DR.dat'
f2 = 'pp_tot_BG.dat'

sigTot(p) = 48.0 + 0.522*log(p)**2 - 4.51*log(p)
sigEl(p) = 11.9 + 26.9*p**(-1.21) + 0.169*log(p)**2 - 1.85*log(p)

################################################################################

set output 'pp_tot_plab.eps'

set xrange [0:8]
set xlabel "p_{lab} [GeV]"
set mxtics 2

plot \
     d1 u 2:5:3:4:($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 7 ps 0.6 t 'data (tot.)', \
     d2 u 2:5:3:4:($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 4 ps 0.6 t 'data (elast.)', \
     f1 u 2:3             lt 1 lw w3 t 'tot', \
     f1 u 2:6             lt 2 lw w t 'NN', \
     f1 u 2:7             lt 3 lc 3 lw w t 'N{/Symbol D}', \
     f1 u 2:4             lt 4 lc 3 lw w t 'NN*', \
     f1 u 2:5             lt 5 lc 3 lw w t 'N{/Symbol D}*', \
     f1 u 2:39            lt 3 lc 4 lw w t '{/Symbol DD}', \
     f1 u 2:37            lt 4 lc 4 lw w t '{/Symbol D}N*', \
     f1 u 2:38            lt 5 lc 4 lw w t '{/Symbol DD}*', \
     f2 u 2:($4+$5+$6+$7) lt 7 lw w t 'NN{/Symbol p}', \
     f2 u 2:($8+$9+$10)   lt 8 lw w t '{/Symbol w/f}', \
     d1 u 2:5:3:4:($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 7 ps 0.6 not, \
     d2 u 2:5:3:4:($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 4 ps 0.6 not, \
     x>4.?sigTot(x):1/0 not

################################################################################

set output 'pp_tot_Ekin.eps'

set xrange [0:7]
set xlabel "E_{kin} [GeV]"

mN = 0.938
Ekin(p)= (p**2+mN**2)**0.5 - mN

plot \
     d1 u (Ekin($2)):5:(Ekin($3)):(Ekin($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 7 ps 0.6 t 'data (tot.)', \
     d2 u (Ekin($2)):5:(Ekin($3)):(Ekin($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 4 ps 0.6 t 'data (elast.)', \
     f1 u (Ekin($2)):3             lt 1 lw w3 t 'tot', \
     f1 u (Ekin($2)):6             lt 2 lw w t 'NN', \
     f1 u (Ekin($2)):7             lt 3 lc 3 lw w t 'N{/Symbol D}', \
     f1 u (Ekin($2)):4             lt 4 lc 3 lw w t 'NN*', \
     f1 u (Ekin($2)):5             lt 5 lc 3 lw w t 'N{/Symbol D}*', \
     f1 u (Ekin($2)):39            lt 3 lc 4 lw w t '{/Symbol DD}', \
     f1 u (Ekin($2)):37            lt 4 lc 4 lw w t '{/Symbol D}N*', \
     f1 u (Ekin($2)):38            lt 5 lc 4 lw w t '{/Symbol DD}*', \
     f2 u (Ekin($2)):($4+$5+$6+$7) lt 7 lw w t 'NN{/Symbol p}', \
     f2 u (Ekin($2)):($8+$9+$10)   lt 8 lw w t '{/Symbol w/f}', \
     f2 u (Ekin($2)):11            lt 9 lw w t 'BYK', \
     d1 u (Ekin($2)):5:(Ekin($3)):(Ekin($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 7 ps 0.6 not, \
     d2 u (Ekin($2)):5:(Ekin($3)):(Ekin($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 4 ps 0.6 not

################################################################################

set output 'pp_tot_srts.eps'

set xrange [1.8:4]
set xlabel "sqrt(s) [GeV]"
set mxtics 5

mN = 0.938
srts(p)= (((p**2+mN**2)**0.5 + mN)**2 - p**2)**0.5
plab(s) = sqrt( ((s**2-2*mN**2)/(2*mN))**2 - mN**2 )

set key samplen 2
#maxrows 6

plot \
     d1 u (srts($2)):5:(srts($3)):(srts($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 7 ps 0.6 t 'data (tot.)', \
     d2 u (srts($2)):5:(srts($3)):(srts($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 4 ps 0.6 t 'data (elast.)', \
     f1 u (srts($2)):3             lt 1 lw w3 t 'tot', \
     f1 u (srts($2)):6             lt 2 lw w t 'NN', \
     f1 u (srts($2)):7             lt 3 lc 3 lw w t 'N{/Symbol D}', \
     f1 u (srts($2)):4             lt 4 lc 3 lw w t 'NN*', \
     f1 u (srts($2)):5             lt 5 lc 3 lw w t 'N{/Symbol D}*', \
     f1 u (srts($2)):39            lt 3 lc 4 lw w t '{/Symbol DD}', \
     f1 u (srts($2)):37            lt 4 lc 4 lw w t '{/Symbol D}N*', \
     f1 u (srts($2)):38            lt 5 lc 4 lw w t '{/Symbol DD}*', \
     f2 u (srts($2)):($4+$5+$6+$7) lt 7 lw w t 'NN{/Symbol p}', \
     f2 u (srts($2)):($8+$9+$10)   lt 8 lw w t '{/Symbol w/f}', \
     f2 u (srts($2)):11            lt 9 lw w t 'BYK', \
     d1 u (srts($2)):5:(srts($3)):(srts($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 7 ps 0.6 not, \
     d2 u (srts($2)):5:(srts($3)):(srts($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 4 ps 0.6 not, \
     x>3.0?sigTot(plab(x)):1/0 not

################################################################################

set output 'pp_inelast.eps'

set xrange [1.9:4]
set yrange [0:40]

set key top left
#maxrows 5

d1 = 'pp_inelast.dat'

sigInel(p) = sigTot(p) - sigEl(p)

plot \
     d1 u (srts($2)):5:(srts($3)):(srts($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 7 ps 0.6 t 'data (inel.)', \
     f1 u (srts($2)):($3-$6)       lt 1 lw w3 t 'inel', \
     f1 u (srts($2)):7             lt 2 lc 2 lw w t 'N{/Symbol D}', \
     f1 u (srts($2)):4             lt 3 lc 2 lw w t 'NN*', \
     f1 u (srts($2)):5             lt 4 lc 2 lw w t 'N{/Symbol D}*', \
     f1 u (srts($2)):39            lt 2 lc 3 lw w t '{/Symbol DD}', \
     f1 u (srts($2)):37            lt 3 lc 3 lw w t '{/Symbol D}N*', \
     f1 u (srts($2)):38            lt 4 lc 3 lw w t '{/Symbol DD}*', \
     f2 u (srts($2)):($4+$5+$6+$7) lt 5 lw w t 'NN{/Symbol p}', \
     f2 u (srts($2)):($8+$9+$10)   lt 7 lw w t '{/Symbol w/f}', \
     f2 u (srts($2)):11            lt 8 lw w t 'BYK', \
     d1 u (srts($2)):5:(srts($3)):(srts($4)):($5+$6):($5-$7) w xyerrorbars lt 1 lc 0 pt 7 ps 0.6 not, \
     x>3.0?sigInel(plab(x)):1/0 not

################################################################################

#!epstopdf pp_tot_plab.eps
#!epstopdf pp_tot_Ekin.eps
#!epstopdf pp_tot_srts.eps