reset

set macros

fOut = "sigmaStar.Q2_01.eps"
Q2=0.100


set xrange [1.0:2.0]
set log y
set yrange [1:1000]

set lmargin 0
set bmargin 0
set tmargin 0
set rmargin 0



set terminal postscript  eps enhanced color 'Helvetica' 24 size 6.0,3.0 dl 2
set output fOut

set multiplot layout 1,2 scale 1.0,1.0 offset 0.5,0

set xlabel "W [GeV]" offset graph 0.5,0
set ylabel "{/Symbol s}^* [{/Symbol m}b]"

set format x
set format y

set xtics 1.0,0.2,1.95
set mxtics 2

set label 1 "{/=72 {/Symbol s}_L}" at graph 0.92, graph 0.15 right
set label 2 "{/=30 Q^2=0.1 GeV^2}" at graph 0.3, graph 0.92 center

plot "fort.222" u 1:($2==Q2?$3:1/0) w lp pt 7 ps 1.5 lc 0 t "Bosted",\
     "fort.122" u 1:($2==Q2?$3:1/0) w l lw 9 lt 1 lc 1 t "total",\
     "fort.122" u 1:($2==Q2?($6>0?$6:0):1/0) w l lw 5 lt 1 lc 2 t "1pi back",\
     "fort.122" u 1:($2==Q2?($6<0?-$6:0):1/0) w l lw 1 lt 2 lc 2 t "",\
     "fort.122" u 1:($2==Q2?($7>0?$7:0):1/0) w l lw 5 lt 1 lc 3 t "2pi back",\
     "fort.122" u 1:($2==Q2?($7<0?-$7:0):1/0) w l lw 1 lt 2 lc 3 t "",\
     "fort.122" u 1:($2==Q2?($8>0?$8:0):1/0) w l lw 5 lt 1 lc 4 t "DIS",\
     "fort.122" u 1:($2==Q2?($8<0?-$8:0):1/0) w l lw 1 lt 2 lc 4 t ""


set xlabel ""
set ylabel ""
set format x
set format y ""

set xtics 1.0,0.2,2.0
set mxtics 2

set label 1 "{/=72 {/Symbol s}_T}" at graph 0.92, graph 0.15 right
unset label 2

plot "fort.222" u 1:($2==Q2?$4:1/0) w lp pt 7 ps 2 lc 0 t "",\
     "fort.122" u 1:($2==Q2?$9:1/0) w l lw 9 lt 1 lc 1 t "",\
     "fort.122" u 1:($2==Q2?($12>0?$12:0):1/0) w l lw 5 lt 1 lc 2 t "",\
     "fort.122" u 1:($2==Q2?($12<0?-$12:0):1/0) w l lw 1 lt 3 lc 2 t "",\
     "fort.122" u 1:($2==Q2?($13>0?$13:0):1/0) w l lw 5 lt 1 lc 3 t "",\
     "fort.122" u 1:($2==Q2?($13<0?-$13:0):1/0) w l lw 1 lt 3 lc 3 t "",\
     "fort.122" u 1:($2==Q2?($14>0?$14:0):1/0) w l lw 5 lt 1 lc 4 t "",\
     "fort.122" u 1:($2==Q2?($14<0?-$14:0):1/0) w l lw 1 lt 3 lc 4 t ""


unset multiplot


set output
set terminal x11

! fixbb @fOut