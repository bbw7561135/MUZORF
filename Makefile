CC = pgcc

ifeq ($(OS),Windows_NT)
RM=del
else
RM=rm
endif

shockrad: shockrad-mfa-disord.c
	$(CC) -g -O2 -mp shockrad-mfa-disord.c -lm -o shockrad-mfa-disord.exe

all: shockrad-mfa-disord

clean:
	$(RM) shockrad-mfa-disord.exe *~ output.txt *.dat

