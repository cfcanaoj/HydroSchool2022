set terminal x11
set xlabel "x"
dtout = 0.005

set yrange [0.0:1.2]
set xrange [-0.5:0.5]

do for [i = 1:40 ] {
  inpfile=sprintf("snap/sod%05d.dat",i)

  # extract time from the first line
  command=sprintf("awk 'NR==1{print($2)}' %s", inpfile)
  time=system(command)

  set title sprintf("time = %s",time) 
  plot inpfile u 1:2 ti "numerical" w lp ,\
       "sod_ana.dat" u ($1*time/0.2):2 ti "exact solution" w l
  pause 0.1
}
