reset

set log zcb

set xrange [:1.0]
set yrange [:1.8]
set zrange [1e-3:]
set cbrange [1e-3:1]



set style data pm3d
set pm3d

splot "fort.141"
#splot "fort.142"
replot "fort.142"
