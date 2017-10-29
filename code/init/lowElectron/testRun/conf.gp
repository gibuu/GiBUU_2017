 reset


#linetyes fuer linien
set style line 1 lc rgb 'red'          lw 3 lt 1 pt 1 ps 1.5
set style line 2 lc rgb 'blue'         lw 3 lt 2 pt 1 ps 1.5
set style line 3 lc rgb 'black'        lw 3 lt 4 pt 1 ps 1.5  #schwarz
set style line 4 lc rgb '#32CD32'      lw 3 lt 5 pt 1 ps 1.5   #grünlich
set style line 5 lc rgb 'dark-orange'  lw 3 lt 7 pt 1 ps 1.5
set style line 6 lc rgb '#FF00FF' lw 2 lt 8 pt 1 ps 1.5      #Fuchsia
set style line 7 lc rgb 'brown'     lw 2 lt 3 pt 1 ps 1.5
set style line 8 lw 2 pt 1 ps 1.5
set style line 9 lw 2 pt 1 ps 1.5

#linetype fuer data points
set style line 10 lc rgb 'black' lt 1 lw 2 pt 7 ps 1.6
set style line 11 lc rgb 'black' lt 1 lw 2 pt 6 ps 1.6
set style line 12 lc rgb 'black' lt 1 lw 2 pt 4 ps 1.6
set style line 13 lc rgb 'black' lt 2 pt 6  lw 2  ps 1.3
set style line 14 lt 1 lw 2 pt 9 ps 1.3
set style line 15 lt 1 lw 2 pt 4 ps 1.3
set style line 16 lt 1 lw 2 pt 3 ps 1.3


set style line 21 lc rgb 'black'     lw 1 lt 1 pt 1 ps 1.5


set key spacing 1.15
set key samplen 3
set xzeroaxis


to='{/Symbol \256}'
pi='{/Symbol p}'
rho='{/Symbol r}'
set term x11
set terminal post eps enhanced color blacktext dashed "Helvetica" 18 lw 2 dl 2 rounded size 15cm, 10 cm
