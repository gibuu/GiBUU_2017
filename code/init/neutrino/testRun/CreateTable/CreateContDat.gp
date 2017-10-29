reset
set log zcb
set cont base
set cntrparam levels discrete 0.2
set table "ContDat.dat"
unset surface
splot "Neutrino.h2D.X_Y.dat" u 1:2:($3*1e11) w l
unset table