reset

set log y
set yrange [1e-5:]

plot "mass2.1.dat" w l lt 1 lw 6 lc 1 t "rho",\
     "mass1.1.dat" w l lt 1 lw 2 lc 0 t "",\
     "mass2.2.dat" w l lt 2 lw 6 lc 2 t "omega",\
     "mass1.2.dat" w l lt 2 lw 2 lc 0 t "",\
     "mass2.3.dat" w l lt 3 lw 6 lc 3 t "phi",\
     "mass1.3.dat" w l lt 3 lw 2 lc 0 t "",\
     "mass2.4.dat" w l lt 4 lw 6 lc 4 t "K*",\
     "mass2.5.dat" w l lt 6 lw 6 lc 6 t "D*s+"