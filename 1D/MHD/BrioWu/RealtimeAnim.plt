set terminal x11
set xlabel "x"
dtout = 0.005

set yrange [0.0:1.2]
set xrange [-0.5:0.5]

do for [i = 1:20 ] {

  inpfile=sprintf("snap/bw2nd_hlld%05d.dat",i)

  # extract time from the first line
  command=sprintf("awk 'NR==1{print($2)}' %s", inpfile)
  time=system(command)

  set title sprintf("time = %s",time) 
  plot inpfile u 1:2 ti "numerical" w lp ,\
       "Brio-Wu-Refsol.dat" u ($2*time/0.1):3 ti "exact solution" w l
  pause 0.1
}
