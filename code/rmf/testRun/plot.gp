
set terminal pdf

set style data lines
f = 'RMF_set9.dat'

##########

set output 'eos.pdf'

plot f u 1:5

##########

set output 'mstar.pdf'

mN = 0.938
rho0 = 0.168
ft(x) = mN * ( 1.-1./(1.+A*(x/rho0)**B) )

fit ft(x) f using 1:(mN-$4) via A, B

plot f u 1:4, mN-ft(x)
