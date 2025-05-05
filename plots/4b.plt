set lmargin 22
set xlabel '|x-y|' font ',22'
set ylabel 'C(x-y)' font ',22' offset -8,0 rotate by 0
set title 'Correlation Function for L=50' font ',22'
set key font ',20'
set xtics font ", 22"
set ytics font ", 22"
set key top center
#set rmargin 45
#set key at screen 1, graph 1

do for [i=11:20] {
set style line i pt i lw 2
}


#A=2
#b=2
#f(x)=A*cosh((x-32)/b)
#ta=2
#tb=1
ta=0.000000001
tb=2

g(x,ta,tb)=ta*cosh((x-25)/tb)
array tmpA[21]
array tmpB[21]
array tmpAerr[21]
array tmpBerr[21]
array chi2[21]
array column1[21]
array column2[21]
array TTT[21]

#fit f(x) 'corrfunc.dat' u 1:2:13 via A, b
#plot 'corrfunc.dat' u 1:2:13 with errorbars, f(x)

do for [i=11:20] {
    TTT[i]=0.1+1.9*(i-1)/19.
    column1[i]=i+1
    column2[i]=i+21
    fit g(x,ta,tb) '../data/corrfunc.dat' using 1:column1[i]:column2[i] via ta, tb
    tmpA[i]=ta
    tmpB[i]=tb
    tmpAerr[i]=ta_err
    tmpBerr[i]=tb_err
    chi2[i]=(FIT_STDFIT*FIT_STDFIT)
}


plot for [i=11:20] '../data/corrfunc.dat' using 1:column1[i]:column2[i] w errorbars linestyle i notitle ,for [i=11:20] g(x,tmpA[i],tmpB[i]) linestyle i title sprintf("T=%.2f, ξ=%.2f±%.2f, χ/dof=%.2f",TTT[i],tmpB[i],tmpBerr[i],chi2[i])

pause -1
