set colors classic

set xlabel "Q^2 [GeV^2]"
set ylabel "FF"

invmass=1.2

plot [0:1] 'deltaFF.dat' u 1:($2==invmass ? $3 : 1/0) w l t 'c3v', '' u 1:($2==invmass ? $4 : 1/0) w l t 'c4v', '' u 1:($2==invmass ? $5 : 1/0) w l t 'c5v',  '' u 1:($2==invmass ? $6 : 1/0) w l t 'c6v', '' u 1:($2==invmass ? $7 : 1/0) w l t 'c3a', '' u 1:($2==invmass ? $8 : 1/0) w l t 'c4a', '' u 1:($2==invmass ? $9 : 1/0) w l t 'c5a', '' u 1:($2==invmass ? $10 : 1/0) w l t 'c6a'
