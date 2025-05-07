set lmargin 20
set xlabel "T" font ',18'
set ylabel "c_4=(3<Q^2>^2-<Q^4>)/L" font ',18' offset -8,0
set xrange [0:2.1]
set yrange [-0.15:0.15]
set key font ',18'
set xtics font ',18'
set ytics font ',18'
plot '../data/c4.dat' w errorbars lw 2 notitle
pause -1
