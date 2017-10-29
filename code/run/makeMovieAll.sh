grep 'All' povray.sh > povrayAll.sh
. povrayAll.sh
convert -delay 20 Movie_All_*.png movie.gif

