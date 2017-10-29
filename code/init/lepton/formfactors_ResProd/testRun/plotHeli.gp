reset

#linetyes fuer linien
set style line 1 lt 1 lw 2 pt 1 ps 1.5
set style line 2 lt 9 lw 2 pt 1 ps 1.5
set style line 3 lt 2 lw 2 pt 1 ps 1.5
set style line 4 lt 4 lw 2 pt 1 ps 1.5
set style line 5 lw 2 pt 1 ps 1.5
set style line 6 lw 2 pt 1 ps 1.5
set style line 7 lw 2 pt 1 ps 1.5
set style line 8 lw 2 pt 1 ps 1.5
set style line 9 lw 2 pt 1 ps 1.5

#linetype fuer data points
set style line 10 lt 1 lw 1.5 pt 7 ps 1.5
set style line 11 lt 9 lw 1.5 pt 6 ps 1.5
set style line 12 lt 9 lw 1.5 pt 5 ps 1.5
set style line 13 lt 1 lw 1.5 pt 2 ps 1.5
set style line 14 lt 7 lw 1.5 pt 9 ps 1.9
set style line 15 lt 1 lw 1.5 pt 4 ps 1.5
set style line 16 lt 1 lw 1.5 pt 3 ps 1.5

set key spacing 1.23
set key samplen	2
set xzeroaxis

set terminal post enhanced eps color blacktext dashed 22 lw 2 dl 2.5


#compare with Lalakulich

set key right bottom
set xlabel 'Q^2 [GeV^2]'
set ylabel 'helicity ampl. [GeV^{-1/2}]'
set xrange [0:4]
set yrange [:]

set output 'helicityAmplitudes_Delta_proton_compareLalakulich.eps'
set title 'Delta proton'
plot "fort.31" u 1:2 w l t 'A32 Lala.', "" u 1:5 w p t 'A32 Maid03',  "" u 1:8 w p t 'A32 Maid05', \
     ""  u 1:3 w l t 'A12 Lala.', "" u 1:6 w p t 'A12 Maid03', "" u 1:9 w p t 'A12 Maid05', \
     ""  u 1:4 w l t 'S12 Lala.', "" u 1:7 w p t 'S12 Maid03', "" u 1:10 w p t 'S12 Maid05'

set output 'helicityAmplitudes_Delta_neutron_compareLalakulich.eps'
set title 'Delta neutron'
plot "fort.30" u 1:2 w l t 'A32 Lala.', "" u 1:5 w p t 'A32 Maid03',  "" u 1:8 w p t 'A32 Maid05', \
     ""  u 1:3 w l t 'A12 Lala.', "" u 1:6 w p t 'A12 Maid03', "" u 1:9 w p t 'A12 Maid05', \
     ""  u 1:4 w l t 'S12 Lala.', "" u 1:7 w p t 'S12 Maid03', "" u 1:10 w p t 'S12 Maid05'


set output 'helicityAmplitudes_P11_proton_compareLalakulich.eps'
set title 'P11 proton'
plot "fort.11" u 1:2 w l t 'A32 Lala.', "" u 1:5 w p t 'A32 Maid03',  "" u 1:8 w p t 'A32 Maid05', \
     ""  u 1:3 w l t 'A12 Lala.', "" u 1:6 w p t 'A12 Maid03', "" u 1:9 w p t 'A12 Maid05', \
     ""  u 1:4 w l t 'S12 Lala.', "" u 1:7 w p t 'S12 Maid03', "" u 1:10 w p t 'S12 Maid05'

set output 'helicityAmplitudes_P11_neutron_compareLalakulich.eps'
set title 'P11 neutron'
plot "fort.10" u 1:2 w l t 'A32 Lala.', "" u 1:5 w p t 'A32 Maid03',  "" u 1:8 w p t 'A32 Maid05', \
     ""  u 1:3 w l t 'A12 Lala.', "" u 1:6 w p t 'A12 Maid03', "" u 1:9 w p t 'A12 Maid05', \
     ""  u 1:4 w l t 'S12 Lala.', "" u 1:7 w p t 'S12 Maid03', "" u 1:10 w p t 'S12 Maid05'


set output 'helicityAmplitudes_S11_proton_compareLalakulich.eps'
set title 'S11 proton'
plot "fort.41" u 1:2 w l t 'A32 Lala.', "" u 1:5 w p t 'A32 Maid03',  "" u 1:8 w p t 'A32 Maid05', \
     ""  u 1:3 w l t 'A12 Lala.', "" u 1:6 w p t 'A12 Maid03', "" u 1:9 w p t 'A12 Maid05', \
     ""  u 1:4 w l t 'S12 Lala.', "" u 1:7 w p t 'S12 Maid03', "" u 1:10 w p t 'S12 Maid05'

set output 'helicityAmplitudes_S11_neutron_compareLalakulich.eps'
set title 'S11 neutron'
plot "fort.40" u 1:2 w l t 'A32 Lala.', "" u 1:5 w p t 'A32 Maid03',  "" u 1:8 w p t 'A32 Maid05', \
     ""  u 1:3 w l t 'A12 Lala.', "" u 1:6 w p t 'A12 Maid03', "" u 1:9 w p t 'A12 Maid05', \
     ""  u 1:4 w l t 'S12 Lala.', "" u 1:7 w p t 'S12 Maid03', "" u 1:10 w p t 'S12 Maid05'

set key right top
set output 'helicityAmplitudes_D13_proton_compareLalakulich.eps'
set title 'D13 proton'
plot "fort.21" u 1:2 w l t 'A32 Lala.', "" u 1:5 w p t 'A32 Maid03',  "" u 1:8 w p t 'A32 Maid05', \
     ""  u 1:3 w l t 'A12 Lala.', "" u 1:6 w p t 'A12 Maid03', "" u 1:9 w p t 'A12 Maid05', \
     ""  u 1:4 w l t 'S12 Lala.', "" u 1:7 w p t 'S12 Maid03', "" u 1:10 w p t 'S12 Maid05'

set key right bottom
set output 'helicityAmplitudes_D13_neutron_compareLalakulich.eps'
set title 'D13 neutron'
plot "fort.20" u 1:2 w l t 'A32 Lala.', "" u 1:5 w p t 'A32 Maid03',  "" u 1:8 w p t 'A32 Maid05', \
     ""  u 1:3 w l t 'A12 Lala.', "" u 1:6 w p t 'A12 Maid03', "" u 1:9 w p t 'A12 Maid05', \
     ""  u 1:4 w l t 'S12 Lala.', "" u 1:7 w p t 'S12 Maid03', "" u 1:10 w p t 'S12 Maid05'



     
#without Lala

set output 'helicityAmplitudes_Delta_proton.eps'
set title 'Delta proton'
plot "fort.31" u 1:5 w l t 'A32 Maid03',  "" u 1:8 w l t 'A32 Maid05', \
     ""  u 1:6 w l t 'A12 Maid03', "" u 1:9 w l t 'A12 Maid05', \
     ""  u 1:7 w l t 'S12 Maid03', "" u 1:10 w l t 'S12 Maid05'

set output 'helicityAmplitudes_Delta_neutron.eps'
set title 'Delta neutron'
plot "fort.30" u 1:5 w l t 'A32 Maid03',  "" u 1:8 w l t 'A32 Maid05', \
     ""  u 1:6 w l t 'A12 Maid03', "" u 1:9 w l t 'A12 Maid05', \
     ""  u 1:7 w l t 'S12 Maid03', "" u 1:10 w l t 'S12 Maid05'


set output 'helicityAmplitudes_P11_proton.eps'
set title 'P11 proton'
plot "fort.11" u 1:5 w l t 'A32 Maid03',  "" u 1:8 w l t 'A32 Maid05', \
     ""  u 1:6 w l t 'A12 Maid03', "" u 1:9 w l t 'A12 Maid05', \
     ""  u 1:7 w l t 'S12 Maid03', "" u 1:10 w l t 'S12 Maid05'

set output 'helicityAmplitudes_P11_neutron.eps'
set title 'P11 neutron'
plot "fort.10" u 1:5 w l t 'A32 Maid03',  "" u 1:8 w l t 'A32 Maid05', \
     ""  u 1:6 w l t 'A12 Maid03', "" u 1:9 w l t 'A12 Maid05', \
     ""  u 1:7 w l t 'S12 Maid03', "" u 1:10 w l t 'S12 Maid05'


set output 'helicityAmplitudes_S11_proton.eps'
set title 'S11 proton'
plot "fort.41" u 1:5 w l t 'A32 Maid03',  "" u 1:8 w l t 'A32 Maid05', \
     ""  u 1:6 w l t 'A12 Maid03', "" u 1:9 w l t 'A12 Maid05', \
     ""  u 1:7 w l t 'S12 Maid03', "" u 1:10 w l t 'S12 Maid05'

set output 'helicityAmplitudes_S11_neutron.eps'
set title 'S11 neutron'
plot "fort.40" u 1:5 w l t 'A32 Maid03',  "" u 1:8 w l t 'A32 Maid05', \
     ""  u 1:6 w l t 'A12 Maid03', "" u 1:9 w l t 'A12 Maid05', \
     ""  u 1:7 w l t 'S12 Maid03', "" u 1:10 w l t 'S12 Maid05'

set key right top
set output 'helicityAmplitudes_D13_proton.eps'
set title 'D13 proton'
plot "fort.21" u 1:5 w l t 'A32 Maid03',  "" u 1:8 w l t 'A32 Maid05', \
     ""  u 1:6 w l t 'A12 Maid03', "" u 1:9 w l t 'A12 Maid05', \
     ""  u 1:7 w l t 'S12 Maid03', "" u 1:10 w l t 'S12 Maid05'

set key right bottom
set output 'helicityAmplitudes_D13_neutron.eps'
set title 'D13 neutron'
plot "fort.20" u 1:5 w l t 'A32 Maid03',  "" u 1:8 w l t 'A32 Maid05', \
     ""  u 1:6 w l t 'A12 Maid03', "" u 1:9 w l t 'A12 Maid05', \
     ""  u 1:7 w l t 'S12 Maid03', "" u 1:10 w l t 'S12 Maid05'


     

!fixbb *.eps
