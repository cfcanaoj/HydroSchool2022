set terminal x11
set xlabel "x"

set yrange [-1.1e-1:1.1e-1]
set xrange [-0.5:0.5]

do for [i = 1:40 ] {

  inpfile=sprintf("lax/snap%05d.dat",i)

  # extract time from the first line
  command=sprintf("awk 'NR==1{print($2)}' %s", inpfile)
  time=system(command)

  set title sprintf("time = %s",time) 
  plot inpfile u 1:8 ti "numerical" w lp ,\
       -1e-1*sin(2.0*pi*(x-time)) ti "exact solution" w l 
  pause 0.1
}
