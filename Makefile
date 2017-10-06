NGmerge: stitch_0_7.c stitch_0_7.h
	gcc -g -Wall -std=gnu99 -fopenmp -O2 -o NGmerge stitch_0_7.c -lz
