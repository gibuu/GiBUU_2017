

# CC muon neutrino total cross sections on neutron

! rm *.eps

reset

load "gnu_config.gp"


set xlabel 'E_{/Symbol n} [GeV]'
set ylabel '{/Symbol s} [10^{-38} cm^2]'

set xrange [0.:3.]

set key right bottom


# QE on neutron
set output 'QE_sigma_CC_muon_neutrino_neutron.eps'
plot 'QE_sigma_CC_muon_neutrino.dat' u 1:3 t 'QE' w l


# delta
set output 'P33_1232_sigma_CC_muon_neutrino_neutron.eps'
set label 'P_{33}(1232)'  at 0.2,.25
plot 'P33_1232_sigma_CC_muon_neutrino.dat' u 1:3 t 'MAID 05' w l ls 1,'P33_1232_sigma_CC_muon_neutrino_Lala__.dat' u 1:3 t 'Lalakulich fit' w l ls 2, 'P33_1232_sigma_CC_muon_neutrino_RuS___.dat' u 1:3 t 'Rein Sehgal' w l ls 3
unset label

set key left top

# P11_1440
set output 'P11_1440_sigma_CC_muon_neutrino_neutron.eps'
set label 'P_{11}(1440)'  at 0.2,0.03
plot 'P11_1440_sigma_CC_muon_neutrino.dat' u 1:3 t 'MAID 05' w l ls 1,'P11_1440_sigma_CC_muon_neutrino_Lala__.dat' u 1:3 t 'Lalakulich fit' w l ls 2, 'P11_1440_sigma_CC_muon_neutrino_RuS___.dat' u 1:3 t 'Rein Sehgal' w l ls 3
unset label

# D13_1520
set output 'D13_1520_sigma_CC_muon_neutrino_neutron.eps'
set label 'D_{13}(1520)'  at 0.2,0.06
plot 'D13_1520_sigma_CC_muon_neutrino.dat' u 1:3 t 'MAID 05' w l ls 1,'D13_1520_sigma_CC_muon_neutrino_Lala__.dat' u 1:3 t 'Lalakulich fit' w l ls 2, 'D13_1520_sigma_CC_muon_neutrino_RuS___.dat' u 1:3 t 'Rein Sehgal' w l ls 3
unset label

# S11_1535
set output 'S11_1535_sigma_CC_muon_neutrino_neutron.eps'
set label 'S_{11}(1535)'  at 0.2,0.06
plot 'S11_1535_sigma_CC_muon_neutrino.dat' u 1:3 t 'MAID 05' w l ls 1,'S11_1535_sigma_CC_muon_neutrino_Lala__.dat' u 1:3 t 'Lalakulich fit' w l ls 2, 'S11_1535_sigma_CC_muon_neutrino_RuS___.dat' u 1:3 t 'Rein Sehgal' w l ls 3
unset label

# resonances around 1.7 GeV
set yrange[0:0.05]
set output 'resonances_around_1700_sigma_CC_muon_neutrino_neutron.eps'
plot 'S31_1620_sigma_CC_muon_neutrino.dat' u 1:3 t 'S_{31}(1620)' w l ls 1, 'S11_1650_sigma_CC_muon_neutrino.dat' u 1:3 t 'S_{11}(1650)' w l ls 2,'D15_1675_sigma_CC_muon_neutrino.dat' u 1:3 t 'D_{15}(1675)' w l ls 3, 'F15_1680_sigma_CC_muon_neutrino.dat' u 1:3 t 'F_{15}(1680)' w l ls 4, 'D33_1700_sigma_CC_muon_neutrino.dat' u 1:3 t 'D_{33}(1700)' w l ls 5, 'P13_1720_sigma_CC_muon_neutrino.dat' u 1:3 t 'P_{13}(1720)' w l ls 6

# resonances around 1.9 GeV
set yrange[0:0.005]
set output 'resonances_around_1900_sigma_CC_muon_neutrino_neutron.eps'
plot 'F35_1905_sigma_CC_muon_neutrino.dat' u 1:3 t 'F_{35}(1905)' w l ls 1, 'P31_1910_sigma_CC_muon_neutrino.dat' u 1:3 t 'P_{31}(1910)' w l ls 2,'F37_1950_sigma_CC_muon_neutrino.dat' u 1:3 t 'F_{37}(1950)' w l ls 3



####################################################################

# CC vector and axial form factors

reset

load "gnu_config.gp"


set xlabel 'Q^2 [GeV^2]'
set ylabel 'FF'


set xrange [0:1]

set key right top

set yrange [-6:6]
set output 'QE_CC_FF.eps'
set label 'nucleon'  at 0.4,4
set ytics  (-6,0,6,3.71, -1.27, 1.00,-117.07)
plot 'QE_CC_FF.dat' u 1:2 t 'F1V' w l ls 1, '' u 1:3 t 'F2V' w l ls 2, '' u 1:4 t'FA' w l ls 3, '' u 1:5 t 'FP' w l ls 4
unset label

set key right top

set yrange [-6:6]
set output 'P33_1232_CC_FF.eps'
set label 'P_{33}(1232)'  at 0.4,4
set ytics  (-6,0,6,2.12, 4.10, -5.19, 1.17, 54.27)
plot 'P33_1232_CC_FF.dat' u 1:2 t 'C3V' w l ls 1, '' u 1:3 t 'C4V' w l ls 2, '' u 1:4 t'C5V' w l ls 3, '' u 1:8 t 'C5A' w l ls 4, '' u 1:9 t 'C6A' w l ls 5
#,'P33_1232_CC_FF_Lala__.dat' u 1:2 t 'C3V Lala' w l, '' u 1:3 t 'C4V Lala' w l, '' u 1:4 t'C5V Lala' w l, '' u 1:8 t 'C5A Lala' w l, '' u 1:9 t 'C6A Lala' w l
unset label

set yrange [-6:6]
set output 'P11_1440_CC_FF.eps'
set label 'P_{11}(1440)'  at 0.4,4
set ytics  (-6,0,6, 4.67, -1.12, -0.52, -61.49)
plot 'P11_1440_CC_FF.dat' u 1:2 t 'F1V' w l ls 1, '' u 1:3 t 'F2V' w l ls 2, '' u 1:6 t'FA' w l ls 3, '' u 1:7 t 'FP' w l ls 4
unset label

set yrange [-4:4]
set output 'D13_1520_CC_FF.eps'
set label 'D_{13}(1520)'  at 0.4,3
set ytics  (-4,0,4, -2.98, 2.84, -1.77, -1.24, -57.46)
plot 'D13_1520_CC_FF.dat' u 1:2 t 'C3V' w l ls 1, '' u 1:3 t 'C4V' w l ls 2, '' u 1:4 t'C5V' w l ls 3, '' u 1:8 t 'C5A' w l ls 4, '' u 1:9 t 'C6A' w l ls 5
unset label

set yrange [0:8]
set output 'S11_1535_CC_FF.eps'
set label 'S_{11}(1535)'  at 0.4,6
set ytics  (0,8,7.45,1.08,0.23, 6.65)
plot 'S11_1535_CC_FF.dat' u 1:2 t 'F1V' w l ls 1, '' u 1:3 t 'F2V' w l ls 2, '' u 1:6 t'FA' w l ls 3, '' u 1:7 t 'FP' w l ls 4
unset label

###########################
## pion prod



reset
load "gnu_config.gp"

set output 'D_vptopplp.eps'
set ylabel '{/Symbol s} [10^{-38} cm^2]'
set xlabel 'E_{/Symbol n} [GeV]'
set key right bottom
set ytics 0.2
set xtics 0.5
#set nokey
set label '{/Symbol n}_{/Symbol m} p {/Symbol \256} {/Symbol m}^- {/Symbol p}^+ p' at 0.4,.9
plot [0.3:2.1] [0:1] 'single_pion_sigma_CC_muon_neutrino.dat' u 1:4 t 'full model, M_A=1.05 GeV' w l ls 2,'single_pion_sigma_CC_muon_neutrino_RuS___.dat' u 1:4 t 'Rein Sehgal' w l ls 3,'/home/leitner/dimpfelmoser/neutrino_calculations/results/vacuum_results/DELTA_sigmatot_numu_CC_proton_Paschos05FF_relBW_manley_nocut.data' u 1:2 t 'only Delta, old FF' w l ls 1, '/home/leitner/dimpfelmoser/neutrino_calculations/results/experimental_data/DELTA_sigmatot_numu_pdeltaplpl/Deuterium_ANL_Barish_1979' u 1:4:2:3:($4-$5):($4+$5) t '' w xyerrorbars ls 10 ,'/home/leitner/dimpfelmoser/neutrino_calculations/results/experimental_data/DELTA_sigmatot_numu_pdeltaplpl/Deuterium_ANL_Radecky_1982' u 1:4:2:3:($4-$5):($4+$5) t '' w xyerrorbars ls 12,'/home/leitner/dimpfelmoser/neutrino_calculations/results/experimental_data/DELTA_sigmatot_numu_pdeltaplpl/Deuterium_BNL_Kitagaki_1986' u 1:4:2:3:($4-$5):($4+$5) t '' w xyerrorbars ls 13, '/home/leitner/dimpfelmoser/neutrino_calculations/results/experimental_data/DELTA_sigmatot_numu_pdeltaplpl/Hydrogen_ANL_Champbell_1973' u 1:4:2:3:($4-$5):($4+$5) t '' w xyerrorbars ls 14


reset
load "gnu_config.gp"

set output 'D_vntoppln.eps'
set ylabel '{/Symbol s} [10^{-38} cm^2]'
set xlabel 'E_{/Symbol n} [GeV]'
set key right bottom
set ytics 0.1
set xtics 0.5
#set nokey
set label '{/Symbol n}_{/Symbol m} n {/Symbol \256} {/Symbol m}^- {/Symbol p}^+ n' at 0.4,.27
plot [0.3:2.1] [0:.3]  'single_pion_sigma_CC_muon_neutrino.dat' u 1:($5+$7*.08) t 'full model, w bg' w l ls 4,'single_pion_sigma_CC_muon_neutrino.dat' u 1:($5) t 'full model, M_A=1.05 GeV' w l ls 2,'single_pion_sigma_CC_muon_neutrino_RuS___.dat' u 1:5 t 'Rein Sehgal' w l ls 3,'/home/leitner/dimpfelmoser/neutrino_calculations/results/vacuum_results/DELTA_sigmatot_numu_CC_proton_Paschos05FF_relBW_manley_nocut.data' u 1:($2/9) t'only Delta, old FF' w l ls 1, '/home/leitner/dimpfelmoser/neutrino_calculations/results/experimental_data/DELTA_sigmatot_numu_nun-nppl/Deuterium_ANL_Barish_1979' u 1:4:2:3:($4-$5):($4+$5) t '' w xyerrorbars ls 10,  '/home/leitner/dimpfelmoser/neutrino_calculations/results/experimental_data/DELTA_sigmatot_numu_nun-nppl/Deuterium_ANL_Radecky_1982' u 1:4:2:3:($4-$5):($4+$5) t '' w xyerrorbars ls 12,'/home/leitner/dimpfelmoser/neutrino_calculations/results/experimental_data/DELTA_sigmatot_numu_nun-nppl/Deuterium_BNL_Kitagaki_1986' u 1:4:2:3:($4-$5):($4+$5) t '' w xyerrorbars ls 13



reset
load "gnu_config.gp"
set output 'D_vntopnullp.eps'
set ylabel '{/Symbol s} [10^{-38} cm^2]'
set xlabel 'E_{/Symbol n} [GeV]'
set key right bottom
set ytics 0.1
set xtics 0.5
#set nokey
set label '{/Symbol n}_{/Symbol m} n {/Symbol \256} {/Symbol m}^- {/Symbol p}^0 p' at 0.4,.35
plot [0.3:2.1] [0:.4]  'single_pion_sigma_CC_muon_neutrino.dat' u 1:($6+$7*0.08) t 'full model, w bg' w l ls 4,'single_pion_sigma_CC_muon_neutrino.dat' u 1:($6) t 'full model, M_A=1.05 GeV' w l ls 2,'single_pion_sigma_CC_muon_neutrino_RuS___.dat' u 1:6 t 'Rein Sehgal' w l ls 3,'/home/leitner/dimpfelmoser/neutrino_calculations/results/vacuum_results/DELTA_sigmatot_numu_CC_proton_Paschos05FF_relBW_manley_nocut.data' u 1:($2*2/9) t 'only Delta, old FF' w l ls 1, '/home/leitner/dimpfelmoser/neutrino_calculations/results/experimental_data/DELTA_sigmatot_numu_nun-ppnull/Deuterium_ANL_Barish_1979' u 1:4:2:3:($4-$5):($4+$5) t '' w xyerrorbars ls 10, '/home/leitner/dimpfelmoser/neutrino_calculations/results/experimental_data/DELTA_sigmatot_numu_nun-ppnull/Deuterium_ANL_Radecky_1982' u 1:4:2:3:($4-$5):($4+$5) t '' w xyerrorbars ls 12, '/home/leitner/dimpfelmoser/neutrino_calculations/results/experimental_data/DELTA_sigmatot_numu_nun-ppnull/Deuterium_BNL_Kitagaki_1986' u 1:4:2:3:($4-$5):($4+$5) t '' w xyerrorbars ls 13













!./fixbb *.eps