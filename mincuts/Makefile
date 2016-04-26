CC = gcc
CFLAGS = -c -std=c99 -Wall
LDFLAGS = -o
DEPENDENCIES = graphCuts.h segmenter.h readENVI.h segment_png_io.h config.h
OBJECTS = main.o graphCuts.o segmenter.o readENVI.o segment_png_io.o config.c

%.o: %.c $(DEPENDENCIES)
	$(CC) -c -o $@ $< $(CFLAGS)

all: build clean

debug: CFLAGS += -g
debug: all

build: $(OBJECTS)
	$(CC) $(LDFLAGS) mincuts $(OBJECTS) -lm -lpng

clean:
	rm *.o
