

diro=EXACTsl
dirs=HLLE1st
dirt=HLLC2nd
dirc=cmp


cc=icc -Wall
fc=ifort -extend-source

t-x-rfig= ${diro}/t-x-r.png  ${dirs}/t-x-r.png  ${dirt}/t-x-r.png
x-r-mov = ${diro}/anix-r.mp4 ${dirs}/anix-r.mp4 ${dirt}/anix-r.mp4 ${dirc}/anix-r.mp4

all: ${t-x-rfig} ${x-r-mov}

${dirc}/anix-r.mp4: MakeAnime_cmp.sh x-r_cmp.plt  ${diro}/t*.dat  ${dirs}/t*.dat ${dirt}/t*.dat
	sh MakeAnime_cmp.sh


##################################################
${dirt}/t-x-r.png: ${dirt}/t-x-r.txd t-x-r.plt
	mv ${dirt} snap
	gnuplot t-x-r.plt
	mv snap ${dirt}

${dirt}/anix-r.mp4: ${dirt}/t*.dat MakeAnime.sh x-r.plt 
	mv ${dirt} snap	
	sh MakeAnime.sh
	mv snap ${dirt}

${dirt}/t-x-r.txd: ${dirt}/t*.dat unifysnap.sh
	mv ${dirt} snap	
	sh unifysnap.sh
	mv snap ${dirt}

${dirt}/t*.dat: HLLC2nd.x
	mkdir snap
	./HLLC2nd.x
	mv snap ${dirt}

HLLC2nd.x: HLLC2nd.cpp
	${cc} $< -o HLLC2nd.x

##################################################

${dirs}/t-x-r.png: ${dirs}/t-x-r.txd t-x-r.plt
	mv ${dirs} snap
	gnuplot t-x-r.plt
	mv snap ${dirs}

${dirs}/anix-r.mp4: ${dirs}/t*.dat MakeAnime.sh x-r.plt 
	mv ${dirs} snap	
	sh MakeAnime.sh
	mv snap ${dirs}

${dirs}/t-x-r.txd: ${dirs}/t*.dat unifysnap.sh
	mv ${dirs} snap	
	sh unifysnap.sh
	mv snap ${dirs}

${dirs}/t*.dat: HLLE1st.x
	mkdir snap
	./HLLE1st.x
	mv snap ${dirs}

HLLE1st.x: HLLE1st.cpp
	${cc} $< -o HLLE1st.x

##################################################

${diro}/t-x-r.png: ${diro}/t-x-r.txd t-x-r.plt
	mv ${diro} snap
	gnuplot t-x-r.plt
	mv snap ${diro}

${diro}/anix-r.mp4: ${diro}/t*.dat MakeAnime.sh x-r.plt 
	mv ${diro} snap	
	sh MakeAnime.sh
	mv snap ${diro}

${diro}/t-x-r.txd: ${diro}/t*.dat unifysnap.sh
	mv ${diro} snap
	sh unifysnap.sh
	mv snap ${diro}

${diro}/t*.dat: riemann.x
	mkdir snap
	./riemann.x
	mv snap ${diro}

riemann.x: evolution.f90
	${fc} $< -o riemann.x


##################################################

clean:
	rm -fr *.x ${diro} ${dirs} ${dirt} snap output

