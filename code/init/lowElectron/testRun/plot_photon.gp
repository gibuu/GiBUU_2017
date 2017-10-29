load 'conf.gp'

set output 'photon_res.eps'

gamma='{/Symbol g}'
set multiplot layout 1,2
set xlabel 'E_
plot "gamma_p_to_R_Xsection.new" u 1:2 w l t 'total:'.gamma.' p -> R' ls 1\
, "gamma_p_to_R_Xsection.new" u 1:3 w l t 'only spin 1/2' ls 2\
, "gamma_p_to_R_Xsection.new" u 1:4 w l t 'only spin 3/2' ls 3 \
, "gamma_p_to_R_Xsection.new" u 1:5 w l t 'only spin 5/2' ls 4\
, "gamma_p_to_R_Xsection.dat" u 1:($2+$3) w l t 'total:'.gamma.' p -> R with bug'


plot "gamma_n_to_R_Xsection.new" u 1:2 w l t 'total:'.gamma.' n -> R' ls 1\
, "gamma_n_to_R_Xsection.new" u 1:3 w l t 'only spin 1/2' ls 2\
, "gamma_n_to_R_Xsection.new" u 1:4 w l t 'only spin 3/2' ls 3\
, "gamma_n_to_R_Xsection.new" u 1:5 w l t 'only spin 5/2' ls 4 \
, "gamma_n_to_R_Xsection.dat" u 1:($2+$3) w l t 'total:'.gamma.' n -> R with bug'

unset multiplot
