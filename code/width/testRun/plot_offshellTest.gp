# reset

# set term x11 0	

# set multiplot layout 3,1
# plot [0:] [0:2] 1.076,1.232,'fort.4714' u 1:2 t 'mom' w l, '' u 1:3 t 'pos' w l, '' u 1:5 t 'width' w l, \
# '' u 1:6 t 'energie' w l, '' u 1:7 t 'mass' w l
# splot "fort.101" u 1:2:3 w l
# splot "fort.101" u 5:6:7 w l
# unset multiplot
# pause -1




reset
set term x11 0 
plot [0:] [0:2] 1.076,1.232,'fort.4714' u 1:2 t 'mom' w l, '' u 1:3 t 'pos' w l, '' u 1:5 t 'width' w l, '' u 1:6 t 'energie' w l, '' u 1:7 t 'mass' w l

set term x11 1
set title'p'
plot [][0.5:] 1.076,'fort.4711' u 1:3 t 'H' w lp, '' u 1:4 t 'offmass' w lp, '' u 1:($5*20.) t 'width' w lp, '' u 1:2 t 'abs4' w lp, '' u 6:($1) t '' w l

set term x11 2
set title'r'
plot [][0.5:]  1.076,'fort.4712' u 1:3 t 'H' w lp, '' u 1:4 t 'offmass' w lp, '' u 1:($5*20.) t 'width' w lp, '' u 1:2 t 'abs4' w lp, '' u 6:($1) t '' w l

set term x11 3
set title'E'
plot [][0.5:] 1.076,'fort.4713' u 1:3 t 'H' w lp, '' u 1:4 t 'offmass' w lp, '' u 1:($5*20.) t 'width' w lp, '' u 1:2 t 'abs4' w lp, '' u 6:($1) t '' w l


pause -1
	