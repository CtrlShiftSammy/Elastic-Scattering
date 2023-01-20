set xrange[0:180]
set yrange[-4.5:4.5]
set title "Ar^{3+} at 16 eV"
set xlabel "Scattering Angle (Degrees)"
set ylabel "Differential Cross Section (Atomic Units)"
set format y "10^{%g}"
p 'Iteration Final/Cross_sections_ar3+_code_static_exchange.txt' u 1:(log10($3)) w l lw 2 dashtype 3 lt rgb "#005A32" title "Rutherford's Formula", \
'Iteration Final/Cross_sections_ar3+_khandker.txt' u 1:(log10($2)) w l lw 2 dashtype 3 lt rgb "magenta" title "Khandker et al", \
'Iteration Final/Cross_sections_ar3+_srigengan.txt' u 1:(log10($2)) w l lw 2 dashtype 2 lt rgb "blue" title "Srigengan et al", \
'Iteration Final/Cross_sections_ar3+_srigengan_exp.txt' u 1:(log10($2)) pointtype 12 lw 2 lt rgb "red" title "Srigengan et al", \
'Iteration Final/Cross_sections_ar3+_code_static_exchange.txt' u 1:(log10($2)) w l lw 2 lt rgb "black" title "Present DCS using static and exchange potentials"