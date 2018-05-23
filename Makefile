all:
	g++ -c -std=c++11 -O3 src/main.cpp -o main.o
	g++ -c -std=c++11 -O3 src/timer.cpp -o timer.o
	g++ -o main main.o timer.o
clean:
	rm *.o main
