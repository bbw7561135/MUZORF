CC = gcc

ifeq ($(OS),Windows_NT)
RM=del
else
RM=rm
endif


telescope: telescope.c
	$(CC) -g telescope.c -o telescope.exe

all: telescope

clean:
	$(RM) telescope.exe 
