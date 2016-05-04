CC = gcc
CFLAGS = -c -std=c99 -Wall
LDFLAGS = -o
DEPENDENCIES = graphCuts.h segmenter.h readENVI.h segmentPngIO.h config.h seeds.h
OBJECTS = main.o graphCuts.o segmenter.o readENVI.o segmentPngIO.o config.c seeds.o

%.o: %.c $(DEPENDENCIES)
	$(CC) -c -o $@ $< $(CFLAGS)

all: build clean

debug: CFLAGS += -g
debug: all

build: $(OBJECTS)
	$(CC) -std=c99 $(LDFLAGS) mincuts $(OBJECTS) -lm -lpng

clean:
	rm *.o
