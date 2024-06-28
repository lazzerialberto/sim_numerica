CC = g++
CFLAGS = -Wall -O3 --std=c++11

es_3_1.exe : es_3_1.o random.o motions.o
	$(CC) random.o motions.o es_3_1.o -o es_3_1.exe
es_3_1.o : es_3_1.cpp motions.h random.h FunzioneBase.h
	$(CC) -c es_3_1.cpp -o es_3_1.o $(CFLAGS)
motions.o : motions.cpp motions.h random.h FunzioneBase.h Integral.h
	$(CC) -c motions.cpp -o motions.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe seed.out