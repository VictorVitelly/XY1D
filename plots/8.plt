set lmargin 20
set xlabel "T" font ',18'
set title "Topological susceptibility" font ',20'
set ylabel "χ_t(T)" font ',18' offset -6,0 rotate by 0
set key font ',18'
set xtics font ',18'
set ytics font ',18'
set xrange [0.05:1.05]
plot '../data/q2.dat' u 1:2:3 w errorbars title "L=50" lw 2, '../data/q2.dat' u 1:4:5 w errorbars title "L=100" lw 2
pause -1
