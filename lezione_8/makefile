CC = g++
CFLAGS = -Wall -O3 --std=c++11

es_8_2.exe : es_8_2.o posizione.o utils.o random.o FunzioneBase.h
	$(CC) es_8_2.o posizione.o utils.o random.o -o es_8_2.exe
es_8_2.o : es_8_2.cpp posizione.h utils.h random.h FunzioneBase.h
	$(CC) -c es_8_2.cpp  -o es_8_2.o $(CFLAGS)
es_8_1.exe : es_8_1.o posizione.o utils.o random.o FunzioneBase.h
	$(CC) es_8_1.o posizione.o utils.o random.o -o es_8_1.exe
es_8_1.o : es_8_1.cpp posizione.h utils.h random.h FunzioneBase.h
	$(CC) -c es_8_1.cpp  -o es_8_1.o $(CFLAGS)
posizione.o : posizione.cpp posizione.h
	$(CC) -c posizione.cpp  -o posizione.o $(CFLAGS)
utils.o : utils.cpp utils.h random.h FunzioneBase.h posizione.h
	$(CC) -c utils.cpp -o utils.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe