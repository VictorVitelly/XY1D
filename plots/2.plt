set terminal qt size 1900, 1000
set lmargin 18
set ylabel offset -2,0
set key font ',16'
set xtics font ", 16"
set ytics font ", 16"
set yrange [-1.0:-0.5]

set multiplot layout 4,1 title "T=0.2, L=50" font ',18'
    set ylabel 'E/L' font ',18' offset -5, 0 rotate by 0
    plot '../data/ThC.dat' u 4:5 w lines title 'Cluster'
    plot '../data/ThM.dat' u 4:5 w lines title 'Metropolis'
    set yrange[-4:4]
    set ylabel 'Q' font ',18' offset -2, 0 rotate by 0
    plot '../data/ThC.dat' u 4:6 w lines title 'Cluster'
    set xlabel 'sweeps' font ',16'
    plot '../data/ThM.dat' u 4:6 w lines title 'Metropolis'

pause -1

unset multiplot
