all: f

f:  main.o newtown_ras.o
	gcc -g -Wall -o f main.o newtown_ras.o

main.o: main.c
	gcc -g -Wall -c -o main.o main.c
	
newtown_ras.o: newtown_ras.s
	nasm -g -f elf64 -w+all -o newtown_ras.o newtown_ras.s


clean:
	rm -f *.o
