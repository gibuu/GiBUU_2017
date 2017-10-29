set terminal postscript enhanced color

#############################################################################

set output "RCS.ps"
set title "Real Compton Scattering: {/Symbol g} p {/Symbol \256} R {/Symbol \256} {/Symbol g} p"

set xlabel "W [GeV]"
set xrange [1.0:2.0]
set mxtics 2

set ylabel "{/Symbol s} [{/Symbol m}b]"

w = 2

f = "RCS.dat"
d = "gg_tot.dat"
plot f using 2:3 with lines title "total" lw w, \
     f using 2:4 with lines title "{/Symbol D}(1232)" lw w, \
     f using 2:9 with lines title "D_{13}(1520)" lw w, \
     f using 2:18 with lines title "F_{15}(1680)" lw w, \
     f using 2:23 with lines title "D_{33}(1700)" lw w lc 7, \
     d using 1:2:3 with errorbars title "data" lw w lt 1 lc 0

#############################################################################

set output "gammaN_R.ps"
set title "{/Symbol g} p {/Symbol \256} R"

set xrange [1.0:2.2]

f = "gammaN_R.dat"
plot f using 2:3 with lines title "total" lw w, \
     f using 2:4 with lines title "{/Symbol D}(1232)" lw w, \
     f using 2:9 with lines title "D_{13}(1520)" lw w, \
     f using 2:18 with lines title "F_{15}(1680)" lw w, \
     f using 2:33 with lines title "F_{37}(1950)" lw w, \
     f using 2:23 with lines title "D_{33}(1700)" lw w lc 7

#############################################################################

set output "BR.ps"
set title "BR(R {/Symbol \256} {/Symbol g} p)"

set logscale y
set ylabel "BR"
set format y "10^{%L}"

f = "BR.dat"
plot f using 2:4 with lines title "{/Symbol D}(1232)" lw w, \
     f using 2:9 with lines title "D_{13}(1520)" lw w, \
     f using 2:18 with lines title "F_{15}(1680)" lw w, \
     f using 2:23 with lines title "D_{33}(1700)" lw w lc 7

#############################################################################
