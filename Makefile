CC = pgcc

ifeq ($(OS),Windows_NT)
RM=del
else
RM=rm
endif

shockrad: shockrad.c
	$(CC) -g -mp shockrad.c -lm -o shockrad.exe

all: shockrad

clean:
	$(RM) shockrad.exe *~ output.txt *.dat

