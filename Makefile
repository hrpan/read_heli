all: main
main: main.cc *.h
	g++ -O2 -Wall -g -o main main.cc -lEG `root-config --cflags` `root-config --libs` -lTreePlayer

clean:
	rm ./*~ ./*.o ./main

