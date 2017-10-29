reset

#set view 90,0
#set view 60,350
set view map

set log zcb

#set zrange [0:0.7]
#set cbrange [0:0.7]

set zrange [0.001:10]
set cbrange [0.05:1]

set xlabel "M [GeV]"
set ylabel "P [GeV]"

set cont both
set cntrparam level discrete 0.08,0.100,0.118,0.130,0.15,0.2,0.25,0.3,0.4
unset surface

splot "DeltaFullWidth.gp1.dat" i 3 u 2:1:4 w l

#set xlabel "M [GeV]"
#set ylabel "dens"
#splot "DeltaFullWidth.gp2.dat" i 1 u 2:3:4 w l