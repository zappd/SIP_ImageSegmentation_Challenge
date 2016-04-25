CC = gcc
CFLAGS = -c -std=c99 -Wall
LDFLAGS = -o
DEPENDENCIES = graphCuts.h segmenter.h readENVI.h
OBJECTS = main.o graphCuts.o segmenter.o readENVI.o

%.o: %.c $(DEPENDENCIES)
	$(CC) -c -o $@ $< $(CFLAGS)

all: build clean

debug: CFLAGS += -g
debug: all

build: $(OBJECTS)
	$(CC) $(LDFLAGS) mincuts $(OBJECTS) -lm

clean:
	rm *.o
