all: f

f:  newtown_ras.o
	gcc -g -Wall -o f newtown_ras.o

	
newtown_ras.o: newtown_ras.s
	nasm -g -f elf64 -w+all -o newtown_ras.o newtown_ras.s


clean:
	rm -f *.o
