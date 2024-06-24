CC = g++
CFLAGS = -Wall -O3 --std=c++11

es_1_3.exe : es_1_3.o random.o
	$(CC) random.o es_1_3.o -o es_1_3.exe
es_1_3.o : es_1_3.cpp 
	$(CC) -c es_1_3.cpp -o es_1_3.o $(CFLAGS)
es_1_2.exe : es_1_2.o random.o
	$(CC) random.o es_1_2.o -o es_1_2.exe
es_1_2.o : es_1_2.cpp 
	$(CC) -c es_1_2.cpp -o es_1_2.o $(CFLAGS)
es_1.exe : es_1.o random.o
	$(CC) random.o es_1.o -o es_1.exe
es_1.o : es_1.cpp 
	$(CC) -c es_1.cpp -o es_1.o $(CFLAGS)
main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe seed.out
