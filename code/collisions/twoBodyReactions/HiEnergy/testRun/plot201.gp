reset
set log y
set yrange [1e-12:100]
set xrange [0:10]

set xlabel "pT [GeV]"
set ylabel "E dsigma/d^3p [mb GeV-2]"

plot "fort.201.TRY3" u 1:3 w l lt 1 t "GiBUU noD", "fort.201.TRY4" u 1:3 w l lt 2 t "GiBUU withD"
replot "fort.201.w0.noDecay" u 1:3 w l lw 4 lt 1 t "QYTHIA noD", "fort.201.w0.withDecay" u 1:3 w l lw 4 lt 2 t "QYTHIA withD"
replot "aaa" u 1:2 w l lw 6 t "Pythia6.225 noD"
replot "bbb" u 1:2 w l lw 6 t "Pythia6.225 withD"