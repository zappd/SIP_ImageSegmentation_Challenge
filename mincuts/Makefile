CC = clang
CFLAGS = -c -std=c11 -Wall
LDFLAGS = -o
DEPENDENCIES = graphCuts.h segmenter.h readENVI.h segmentPngIO.h config.h seeds.h
OBJECTS = main.o graphCuts.o segmenter.o readENVI.o segmentPngIO.o config.c seeds.c

%.o: %.c $(DEPENDENCIES)
	$(CC) -c -o $@ $< $(CFLAGS)

all: build clean

debug: CFLAGS += -g
debug: all

build: $(OBJECTS)
	$(CC) $(LDFLAGS) mincuts $(OBJECTS) -lm -lpng

clean:
	rm *.o
