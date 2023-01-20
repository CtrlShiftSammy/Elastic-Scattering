set xrange[0:180]
set yrange[-4.5:4.5]
set title "Ar^{+} at 16 eV"
set xlabel "Scattering Angle (Degrees)"
set ylabel "Differential Cross Section (Atomic Units)"
set format y "10^{%g}"
p 'Iteration Final/Cross_sections_ar+_code_static_exchange.txt' u 1:(log10($3)) w l lw 2 dashtype 3 lt rgb "#005A32" title "Rutherford's Formula", \
'Iteration Final/Cross_sections_ar+_khandker.txt' u 1:(log10($2)) w l lw 2 dashtype 3 lt rgb "magenta" title "Khandker et al", \
'Iteration Final/Cross_sections_ar+_brotton.txt' u 1:(log10($2)) w l lw 2 dashtype 2 lt rgb "blue" title "Brotton et al", \
'Iteration Final/Cross_sections_ar+_brotton_exp.txt' u 1:(log10($2)) pointtype 12 lw 2 lt rgb "red" title "Brotton et al", \
'Iteration Final/Cross_sections_ar+_code_static_exchange.txt' u 1:(log10($2)) w l lw 2 lt rgb "black" title "Present DCS using static and exchange potentials"