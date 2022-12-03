reset

pngflag=1

if(pngflag==1) outflag=1
if(outflag==1) set terminal push

if(pngflag ==1) set terminal pngcairo enhanced font "Helvetica, 18" crop

input="snap/t-x-r.txd"

set pm3d map

set xlabel "x" offset 0.0, 1.0
set xtic offset 0.0,0.5
set xtic 0.2
set mxtics 2

set ylabel "t" offset 1,0
set ytic offset 0.5,0.0
set ytic 0.05
set mytic 5
set yrange [*:0.2]

##########################################
# Space Time Diagram
##########################################
if(pngflag ==1)set output "snap/t-x-r.png"

set palette defined ( 0 "black", 0.1 "blue", 0.5 "green", 1.0 "dark-green" )
set cbrange [*:*]
set title "{/Symbol r}"
set title offset -12.0,-0.8
splot input  u 2:1:3 notitle \

if(pngflag ==1)set output "snap/t-x-p.png"

set palette defined ( 0 "black", 0.1 "blue", 0.2 "red",0.5 "orange", 1.0 "yellow"  )
set cbrange [*:*]
set title "p"

splot input  u 2:1:4 notitle \

if(pngflag ==1)set output "snap/t-x-vr.png"

set palette defined ( 0 "white", 1.0 "red" )
set cbrange [*:*]
set title "v"

splot input  u 2:1:5 notitle \

if(pngflag ==1)set output "snap/t-x-s.png"
set title "s"

set palette defined ( 0 "black", 0.1 "purple",0.2 "blue", 0.8 "red",0.9 "orange", 1.0 "yellow"  )
set cbrange [*:*]

rho0=1.0
p0=1.0
s0=1.0
gamma=1.4
entlopy(rho,p)=1.0/(gamma-1.0)*log(p/p0*(rho0/rho)**(gamma))+s0

splot input  u 2:1:(entlopy($3,$4)) notitle \

##########################################
# Finalize
##########################################

reset
if(outflag==1) set terminal pop
