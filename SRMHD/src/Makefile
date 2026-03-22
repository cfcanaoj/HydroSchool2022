FC      = gfortran
FFLAGS  = -O2
#FFLAGS = -fopenmp

MODS = params_mod.o riemann_mod.o prim_mod.o flux_mod.o

TARGET_2D = srmhd_khi.x
TARGET_1D = main1d.x

.PHONY: all 2d 1d clean debug

all: 2d 1d
2d: $(TARGET_2D)
1d: $(TARGET_1D)

$(TARGET_2D): $(MODS) main2d.o
	$(FC) $(FFLAGS) -o $@ $^

$(TARGET_1D): params_mod.o riemann_mod.o prim_mod.o flux_mod.o main1d.o
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

riemann_mod.o: params_mod.o
prim_mod.o:    params_mod.o
flux_mod.o:    params_mod.o riemann_mod.o prim_mod.o
main2d.o:      $(MODS)
main1d.o:      params_mod.o riemann_mod.o prim_mod.o flux_mod.o

debug: FFLAGS = -O0 -g -fopenmp -Wall -Wextra -fcheck=all -fbacktrace
debug: $(TARGET_2D) $(TARGET_1D)

clean:
	rm -f $(TARGET_2D) $(TARGET_1D) *.o *.mod
