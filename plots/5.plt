set lmargin 15
set xlabel "T" font ',18'
set ylabel "M^2" font ',18' offset -2,0 rotate by 0
set key font ',18'
set xtics font ',18'
set ytics font ',18'
plot '../data/magnet2.dat' u 1:2:3 w errorbars title "L=50" lw 2, '../data/magnet2.dat' u 1:4:5 w errorbars title "L=100" lw 2
pause -1
