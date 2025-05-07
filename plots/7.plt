set lmargin 15
set xlabel "T" font ',18'
set ylabel "ξ(T)" font ',18' offset -4,0 rotate by 0
set title "Correlation length" font ',20'
set key font ',18'
set xtics font ',18'
set ytics font ',18'
set xrange [0.05:1.05]
plot '../data/xi.dat' u 1:2:3 w errorbars title "L=50" lw 2, '../data/xi.dat' u 1:4:5 w errorbars title "L=100" lw 2
pause -1
