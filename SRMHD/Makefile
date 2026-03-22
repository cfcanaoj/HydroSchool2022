.PHONY: all 2d 1d clean debug run2d run1d

all: 2d 1d

2d:
	$(MAKE) -C src 2d

1d:
	$(MAKE) -C src 1d

clean:
	$(MAKE) -C src clean

debug:
	$(MAKE) -C src debug

# Run binaries from project root so that output dirs (ct/, bw1d/) are created here
run2d: 2d
	./src/srmhd_khi.x

run1d: 1d
	./src/main1d.x
