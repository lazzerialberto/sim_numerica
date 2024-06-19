CC = g++
CFLAGS = -Wall -O3 --std=c++11

es_2_2.exe : es_2_2.o posizione.o
	$(CC) es_2_2.o posizione.o -o es_2_2.exe
es_2_2.o : es_2_2.cpp posizione.h
	$(CC) -c es_2_2.cpp -o es_2_2.o $(CFLAGS)
generate_rws_c.exe : generate_rws_c.o posizione.o random.o
	$(CC) generate_rws_c.o posizione.o random.o -o generate_rws_c.exe
generate_rws_c.o : generate_rws_c.cpp posizione.h
	$(CC) -c generate_rws_c.cpp  $(CFLAGS)
generate_rws.exe : generate_rws.o posizione.o random.o
	$(CC) generate_rws.o posizione.o random.o -o generate_rws.exe
generate_rws.o : generate_rws.cpp posizione.h
	$(CC) -c generate_rws.cpp  $(CFLAGS)
posizione.o : posizione.cpp posizione.h
	$(CC) -c posizione.cpp  $(CFLAGS)
es_2_1.exe : es_2_1.o random.o
	$(CC) random.o es_2_1.o -o es_2_1.exe
es_2_1.o : es_2_1.cpp 
	$(CC) -c es_2_1.cpp -o es_2_1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out