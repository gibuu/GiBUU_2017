
################### general stuff ###################

set terminal postscript enhanced color

set xlabel 'E_{/Symbol g} [GeV]'
set mxtics 2

set ylabel '{/Symbol s}_{{/Symbol g}N} [{/Symbol m}b]'
set mytics 2

set style data lines
w = 3

f = 'vecmes001.dat'
f1 = 'vecmes001.dat'
f2 = 'vecmes002.dat'
f3 = 'vecmes003.dat'

###########################################

set output 'gammaN_all.ps'
set xrange [0:10]
#set yrange [0:30]
#set title '{/Symbol g} N {/Symbol \256 r} X'

plot f1 using 1:3        lw w lt 1 lc 1 title '{/Symbol r} N (Effe)', \
     f1 using 1:4        lw w lt 1 lc 2 title '{/Symbol w} N (Effe)', \
     f1 using 1:5        lw w lt 1 lc 3 title '{/Symbol f} N (Effe)', \
     f2 using 1:3        lw w lt 2 lc 1 title '{/Symbol r} N (Pythia)', \
     f2 using 1:4        lw w lt 2 lc 2 title '{/Symbol w} N (Pythia)', \
     f2 using 1:5        lw w lt 2 lc 3 title '{/Symbol f} N (Pythia)', \
     f3 using 1:3        lw w lt 3 lc 1 title '{/Symbol r} N (Donnachie)', \
     f3 using 1:4        lw w lt 3 lc 2 title '{/Symbol w} N (Donnachie)', \
     f3 using 1:5        lw w lt 3 lc 3 title '{/Symbol f} N (Donnachie)'

################### RHO ###################

set output 'gammaN_rhoX.ps'
set xrange [0.5:3.0]
set yrange [0:30]
set title '{/Symbol g} N {/Symbol \256 r} X'

plot f using 1:3        lw w lt 1 lc 2 title '{/Symbol r} N', \
     f using 1:15       lw w lt 2 lc 2 title 'R{/Symbol \256 r} N', \
     f using 1:11       lw w lt 1 lc 3 title '{/Symbol r D}', \
     f using 1:18       lw w lt 2 lc 3 title 'R{/Symbol \256 r D}', \
     f using 1:($3+$11) lw w lt 1 lc 1 title 'total'

################### OMEGA ###################

set output 'gammaN_omegaX.ps'
set xrange [1.0:3.0]
set yrange [0:10]
set title '{/Symbol g} N {/Symbol \256 w} X'

plot f using 1:4        lw w lt 1 lc 2 title '{/Symbol w} N', \
     f using 1:16       lw w lt 2 lc 2 title 'R{/Symbol \256 w} N', \
     f using 1:12       lw w lt 1 lc 3 title '{/Symbol w D}', \
     f using 1:19       lw w lt 2 lc 3 title 'R{/Symbol \256 w D}', \
     f using 1:($4+$12) lw w lt 1 lc 1 title 'total', \
     '../../../../../../buuinput/gammaN_omegaN_saphir.dat' with points lc 0 title 'data'


################### PHI ###################

set output 'gammaN_phiX.ps'
set xrange [1.5:4.0]
set yrange [0:1]
set title '{/Symbol g} N {/Symbol \256 f} X'

plot f using 1:5        lw w lt 1 lc 2 title '{/Symbol f} N', \
     f using 1:17       lw w lt 2 lc 2 title 'R{/Symbol \256 f} N', \
     f using 1:14       lw w lt 1 lc 3 title '{/Symbol f D}', \
     f using 1:20       lw w lt 2 lc 3 title 'R{/Symbol \256 f D}', \
     f using 1:($5+$14) lw w lt 1 lc 1 title 'total'

################### ETA ###################

set output 'gammaN_etaX.ps'
set xrange [0.5:2.0]
set yrange [0:20]
set title '{/Symbol g} N {/Symbol \256 h} X'

plot f using 1:21       lw w lt 2 lc 2 title 'R{/Symbol \256 h} N'


################### total VMD XS ###################

set output 'gammaN_VMD.ps'
set xrange [1.0:3.0]
set yrange [0:150]
set title '{/Symbol g} N {/Symbol \256} X'

plot f using 1:22       lw w lt 2 lc 2 title 'VMD'