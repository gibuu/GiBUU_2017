#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.0 patchlevel 0
#    	last modified Thu Apr 15 14:44:22 CEST 2004
#    	System: Linux 2.6.18.8-0.1-default
#    
#    	Copyright (C) 1986 - 1993, 1998, 2004
#    	Thomas Williams, Colin Kelley and many others
#    
#    	This is gnuplot version 4.0.  Please refer to the documentation
#    	for command syntax changes.  The old syntax will be accepted
#    	throughout the 4.0 series, but all save files use the new syntax.
#    
#    	Type `help` to access the on-line reference manual.
#    	The gnuplot FAQ is available from
#    		http://www.gnuplot.info/faq/
#    
#    	Send comments and requests for help to
#    		<gnuplot-info@lists.sourceforge.net>
#    	Send bugs, suggestions and mods to
#    		<gnuplot-bugs@lists.sourceforge.net>
#    
# set terminal x11 
# set output
unset clip points
set clip one
unset clip two
set bar 1.000000
set border 31 lt -1 lw 1.000
set xdata
set ydata
set zdata
set x2data
set y2data
set timefmt x "%d/%m/%y,%H:%M"
set timefmt y "%d/%m/%y,%H:%M"
set timefmt z "%d/%m/%y,%H:%M"
set timefmt x2 "%d/%m/%y,%H:%M"
set timefmt y2 "%d/%m/%y,%H:%M"
set timefmt cb "%d/%m/%y,%H:%M"
set boxwidth
set style fill empty border
set dummy x,y
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set angles radians
unset grid
set key title ""
set key right top Right noreverse enhanced box linetype -2 linewidth 1.000 samplen 4 spacing 1 width 0 height 0 autotitles
unset label
unset arrow
unset style line
unset style arrow
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis lt -2 lw 1.000
set yzeroaxis lt -2 lw 1.000
set x2zeroaxis lt -2 lw 1.000
set y2zeroaxis lt -2 lw 1.000
set tics in
set ticslevel 0.5
set ticscale 1 0.5
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics border mirror norotate autofreq 
set ytics border mirror norotate autofreq 
set ztics border nomirror norotate autofreq 
set nox2tics
set noy2tics
set cbtics border mirror norotate autofreq 
set title "" 0.000000,0.000000  font ""
set timestamp "" bottom norotate 0.000000,0.000000  ""
set rrange [ * : * ] noreverse nowriteback  # (currently [0.00000:10.0000] )
set trange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set urange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set vrange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set xlabel "" 0.000000,0.000000  font ""
set x2label "" 0.000000,0.000000  font ""
set xrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set x2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set ylabel "" 0.000000,0.000000  font ""
set y2label "" 0.000000,0.000000  font ""
set yrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set y2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set zlabel "" 0.000000,0.000000  font ""
set zrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set cblabel "" 0.000000,0.000000  font ""
set cbrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set zero 1e-08
set lmargin -1
set bmargin -1
set rmargin -1
set tmargin -1
set locale "C"
set pm3d scansautomatic flush begin noftriangles nohidden3d implicit corners2color mean
unset pm3d
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin 0.9,0.2 size 0.1,0.63 bdefault
set loadpath 
set fontpath 
set fit noerrorvariables
MOUSE_X = 0.878043418521572
MOUSE_Y = 0.104575197889182
MOUSE_X2 = 0.836702390766694
MOUSE_Y2 = 0.102408997361478
MOUSE_BUTTON = 1
MOUSE_SHIFT = 0
MOUSE_ALT = 0
MOUSE_CTRL = 0
set xlabel 'E_final'
set ylabel 'dsigma/dOmega/dE [nb/MeV]'
plot "fort.12" u ($1):($2*1000000) w l t 'P_33' 
replot "fort.13" u ($1):($2*1000000) w l t 'P_11(1440)'
replot "fort.14" u ($1):($2*1000000) w l t 'S_11(1535)'
replot "fort.15" u ($1):($2*1000000) w l t 'S_11(1650)'
#replot "fort.16" u ($1):($2*1000000) w l t 'S_11(2090)'
replot "fort.17" u ($1):($2*1000000) w l t 'D_13(1520)'
#replot "fort.18" u ($1):($2*1000000) w l t 'D_13(1700)'
#replot "fort.19" u ($1):($2*1000000) w l t 'D_13(2080)'
replot "fort.20" u ($1):($2*1000000) w l t 'D_15(1675)'
#replot "fort.21" u ($1):($2*1000000) w l t 'G_17(2190)'
#replot "fort.22" u ($1):($2*1000000) w l t 'P_11(1710)'
#replot "fort.23" u ($1):($2*1000000) w l t 'P_11(2100)'
replot "fort.24" u ($1):($2*1000000) w lp t 'P_13(1720)'
#replot "fort.25" u ($1):($2*1000000) w l t 'P_13'
replot "fort.26" u ($1):($2*1000000) w lp t 'F_15(1680)'
#replot "fort.27" u ($1):($2*1000000) w l t 'F_15(2000)'
#replot "fort.28" u ($1):($2*1000000) w l t 'F_17(1990)
replot "fort.29" u ($1):($2*1000000) w lp t 'S31_1620'
#replot "fort.30" u ($1):($2*1000000) w l t  'S31_1900'
replot "fort.31" u ($1):($2*1000000) w lp t 'D33_1700'
#replot "fort.32" u ($1):($2*1000000) w l t 'D33_1940'
#replot "fort.33" u ($1):($2*1000000) w l t 'D35_1930'
#replot "fort.34" u ($1):($2*1000000) w l t  'D35_2350'
#replot "fort.35" u ($1):($2*1000000) w l t  'P31'
replot "fort.36" u ($1):($2*1000000) w lp t  'P31_1910=26'
#replot "fort.37" u ($1):($2*1000000) w l t  'P33_1600=27'
#replot "fort.38" u ($1):($2*1000000) w l t  'P33_1920'
#replot "fort.39" u ($1):($2*1000000) w l t  'F35'
replot "fort.40" u ($1):($2*1000000) w lp t  'F35_1905'
replot "fort.41" u ($1):($2*1000000) w lp t 'F37_1950=31'

replot "fort.100" u ($1):($2*1000000) w lp lt -1 t 'total'

#pause -1

#replot   "~/Documents/DrArbeit/ElectronInduced/paper_Experimente/PRL53_17.jpg.dat.16" u (0.73-$1):($2/16) w p ps 2 pt 6 t 'data'


pause -1

#  !Nucleon Resonances
  #, parameter :: nucleon=1
  #, parameter :: delta=2
  #, parameter :: P11_1440=3
  #, parameter :: S11_1535=4
  #, parameter :: S11_1650=5
  #, parameter :: S11_2090=6
  #, parameter :: D13_1520=7
  #, parameter :: D13_1700=8
  #, parameter :: D13_2080=9
  #, parameter :: D15_1675=10
  #, parameter :: G17_2190=11
  #, parameter :: P11_1710=12
  #, parameter :: P11_2100=13
  #, parameter :: P13_1720=14
  #, parameter :: P13=15
  #, parameter :: F15_1680=16
  #, parameter :: F15_2000=17
  #, parameter :: F17_1990=18
  #, parameter :: S31_1620=19
  #, parameter :: S31_1900=20
  #, parameter :: D33_1700=21
  #, parameter :: D33_1940=22
  #, parameter :: D35_1930=23
  #, parameter :: D35_2350=24
  #, parameter :: P31=25
  #, parameter :: P31_1910=26
  #, parameter :: P33_1600=27
  #, parameter :: P33_1920=28
  #, parameter :: F35=29
  #, parameter :: F35_1905=30
  #, parameter :: F37_1950=31
