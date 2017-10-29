reset

set log y
set xrange [2:4]     # for sqrt(s)
set yrange [0.1:100]
set xlabel "sqrt(s) [GeV]"
set ylabel "{/Symbol s}_{pp} [mb]"
set mxtics 5
set key top left
set key box

set terminal postscript enhanced landscape color 'Helvetica' 18
set output 'Plot23.ps'

w = 5

sqrts0 = 2.658

f(x) = 2.5*((x/sqrts0)**2-1.)**1.47*((x/sqrts0)**2)**-1.11

# Plotting as function of sqrt(s):
plot "fort.23" u 2:3 w l lw w lt 1 t '{/Symbol s}(tot.)', \
     "fort.23" u 2:4 w l lw w lt 2 t '{/Symbol s}(elast.)', \
     "fort.23" u 2:5 w l lw w lt 3 t '{/Symbol s}(wX)', \
     "fort.23" u 2:6 w l lw w lt 4 t '{/Symbol s}(wNN)', \
      f(x) t 'Sibirtsev'
