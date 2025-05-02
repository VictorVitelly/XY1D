set terminal qt size 1300, 900
set lmargin 16
set ylabel 'Freq.' font ',16' rotate by 0
set key font ',16'
binwidth=1.
bin(x,width)=width*floor(x/width)
set boxwidth binwidth
sigma=1
xmean=0.
gauss(x)=exp(-(x-mean)**2 /(2*sigma**2)  )/(sqrt(2*pi)*sigma)
gauss2(x)=exp(-(x-mean2)**2 /(2*sigma2**2)  )/(sqrt(2*pi)*sigma2)
fit gauss(x) '../data/histogram.dat' u 1:2:3 via sigma, mean
chi2=(FIT_STDFIT*FIT_STDFIT)
fit gauss2(x) '../data/histogram.dat' u 4:5:6 via sigma2, mean2
chi22=(FIT_STDFIT*FIT_STDFIT)
#plot 'data/data.dat' using (bin($2,binwidth)):(1.0) smooth freq with boxes notitle, gauss(x)
set multiplot layout 2,1 title "Topological Charge Histogram" font ',18'
    plot '../data/histogram.dat' u 1:2:3 notitle w errorbars, gauss(x) lw 2 title sprintf("<Q>=%.2f±%.2f, σ=%.2f±%.2f, χ^2/dof=%.2f",mean,mean_err,sigma,sigma_err,chi2), '../data/histogram.dat' u 1:2:3 smooth freq with boxes title "T=1.0" lt rgb "purple" lw 2
    set xlabel 'Q' font ',16'
    plot '../data/histogram.dat' u 4:5:6 notitle w errorbars, gauss2(x) lw 2 title sprintf("<Q>=%.3f±%.3f, σ=%.3f±%.3f, χ^2/dof=%.2f",mean2,mean2_err,sigma2,sigma2_err,chi22), '../data/histogram.dat' u 4:5:6 smooth freq with boxes title "T=0.5" lt rgb "purple" lw 2

pause -1
unset multiplot
