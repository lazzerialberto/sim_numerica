CC = g++
CFLAGS = -Wall -O3 --std=c++11

tsp.exe : tsp.o posizione.o utils.o random.o
	$(CC) tsp.o posizione.o utils.o random.o -o tsp.exe
tsp.o : tsp.cpp posizione.h utils.h random.h
	$(CC) -c tsp.cpp  -o tsp.o $(CFLAGS)
posizione.o : posizione.cpp posizione.h
	$(CC) -c posizione.cpp  -o posizione.o $(CFLAGS)
utils.o : utils.cpp utils.h random.h posizione.h
	$(CC) -c utils.cpp -o utils.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe