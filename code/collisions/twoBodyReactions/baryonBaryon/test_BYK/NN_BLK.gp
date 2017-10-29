set terminal pdf enhanced dashed

set output 'NN_BLK.pdf'

set title 'NN -> {/Symbol L}X'

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
      f1 u 1:($4+$8+$12) lt 1 lc 1 lw 3 t 'pp -> {/Symbol L}X', \
      f1 u 1:4           lt 2 lc 1      t 'p{/Symbol L}K^+', \
      f1 u 1:8           lt 3 lc 1      t '{/Symbol D}^{+}{/Symbol L}K^+', \
      f1 u 1:12          lt 4 lc 1      t '{/Symbol D}^{++}{/Symbol L}K^0', \
      f2 u 1:($4+$7+$10+$13) lt 1 lc 2 lw 3 t 'pn -> {/Symbol L}X', \
      f2 u 1:4               lt 2 lc 2      t 'n{/Symbol L}K^+', \
      f2 u 1:10              lt 4 lc 2      t 'p{/Symbol L}K^0', \
      f2 u 1:7               lt 3 lc 2      t '{/Symbol D}^{0}{/Symbol L}K^+', \
      f2 u 1:13              lt 5 lc 2      t '{/Symbol D}^{+}{/Symbol L}K^0'
