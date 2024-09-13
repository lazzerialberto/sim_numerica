CC = g++
CFLAGS = -Wall -O3 --std=c++11

es_5_1.exe : es_5_1.o posizione.o utils.o random.o FunzioneBase.h
	$(CC) es_5_1.o posizione.o utils.o random.o -o es_5_1.exe
es_5_1.o : es_5_1.cpp posizione.h utils.h random.h FunzioneBase.h
	$(CC) -c es_5_1.cpp  -o es_5_1.o $(CFLAGS)
posizione.o : posizione.cpp posizione.h
	$(CC) -c posizione.cpp  -o posizione.o $(CFLAGS)
utils.o : utils.cpp utils.h random.h FunzioneBase.h posizione.h
	$(CC) -c utils.cpp -o utils.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe