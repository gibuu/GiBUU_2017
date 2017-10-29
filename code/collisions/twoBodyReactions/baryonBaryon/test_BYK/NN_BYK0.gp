set terminal pdf enhanced dashed

set output 'NN_BYK0.pdf'

set title 'NN -> K^0X'

set key top left samplen 2

set style data lines

set xlabel "sqrt(s) [GeV]"
set xrange [2.3:3.6]
set mxtics 2

set ylabel "{/Symbol s} [mb]"
set log y
set yrange [2E-4:2]
set mytics 10
set format y "10^{%L}"

f1 = 'N+N+_BarHypKaon_XS.dat'
f2 = 'N+N0_BarHypKaon_XS.dat'

plot \
      f1 u 1:($11+$12+$13+$14) lt 1 lc 1 lw 3 t 'pp -> K^0X', \
      f1 u 1:11                lt 3 lc 1      t 'p{/Symbol S}^+K^0', \
      f1 u 1:12                lt 4 lc 1      t '{/Symbol D}^{++}{/Symbol L}K^0', \
      f1 u 1:13                lt 5 lc 1      t '{/Symbol D}^{++}{/Symbol S}^0K^0', \
      f1 u 1:14                lt 6 lc 1      t '{/Symbol D}^{+}{/Symbol S}^+K^0', \
      f2 u 1:($10+$11+$12+$13+$14+$15) lt 1 lc 2 lw 3 t 'pn -> K^0X', \
      f2 u 1:10                        lt 2 lc 2      t 'p{/Symbol L}^0K^0', \
      f2 u 1:11                        lt 3 lc 2      t 'p{/Symbol S}^0K^0', \
      f2 u 1:12                        lt 4 lc 2      t 'n{/Symbol S}^+K^0', \
      f2 u 1:13                        lt 5 lc 2      t '{/Symbol D}^{+}{/Symbol L}^0K^0', \
      f2 u 1:14                        lt 6 lc 2      t '{/Symbol D}^{0}{/Symbol S}^+K^0', \
      f2 u 1:15                        lt 7 lc 2      t '{/Symbol D}^{+}{/Symbol S}^0K^0'
