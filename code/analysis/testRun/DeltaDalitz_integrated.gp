set terminal postscript eps enhanced dashed color

set output 'Gamma_DeltaDalitz.eps'

set xlabel 'm_{/Symbol D} [GeV]'
set ylabel '{/Symbol G} [GeV]'

set mxtics 2

set log y
set format y "10^{%L}"
set mytics 10
set yrange [1e-8:1e-1]

set key bottom right

set arrow from 1.232,1E-8 to 1.232,1E-1 lt 0 nohead

w = 3

set style data lines

f11 = 'Gamma_DeltaDalitz.11.dat'  # Wolf
f21 = 'Gamma_DeltaDalitz.21.dat'  # Krivo
f25 = 'Gamma_DeltaDalitz.25.dat'
f26 = 'Gamma_DeltaDalitz.26.dat'
f41 = 'Gamma_DeltaDalitz.41.dat'  # Ernst

plot \
     f21 u 1:4 lt 1 lw w t 'real photon (Krivo)', \
     f11 u 1:4 lt 2 lw w t 'real photon (Wolf)', \
     f41 u 1:4 lt 3 lw w t 'real photon (Ernst)', \
     f21 u 1:3 lt 2 lc 1 lw w t 'dilepton (Krivo)', \
     f25 u 1:3 lt 3 lc 1 lw w t 'dilepton (Krivo + Iachello)', \
     f26 u 1:3 lt 4 lc 1 lw w t 'dilepton (Krivo + Ramalho)', \
     f11 u 1:3 lt 3 lc 2 lw w t 'dilepton (Wolf)', \
     f41 u 1:3 lt 4 lc 3 lw w t 'dilepton (Ernst)'
