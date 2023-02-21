set terminal x11
set xlabel "x"

set yrange [0.0:1.2]
set xrange [-0.5:0.5]

do for [i = 1:40 ] {
  inpfile=sprintf("lax/snap%05d.dat",i)

  # extract time from the first line
  command=sprintf("awk 'NR==1{print($2)}' %s", inpfile)
  time=system(command)

  set title sprintf("time = %s",time) 
  plot inpfile u 1:2 ti "numerical" w lp ,\
      "briowu_nonregsol.dat" u (($1-0.5)*time/0.1):2 ti "exact solution" w l
  pause 0.1
}
