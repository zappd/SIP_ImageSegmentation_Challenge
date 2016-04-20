CC = gcc
CFLAGS = -c -Wall
LDFLAGS = -o
DEPENDENCIES = graphCuts.h config.h
OBJECTS = main.o graphCuts.o config.o

%.o: %.c $(DEPENDENCIES)
	$(CC) -c -o $@ $< $(CFLAGS)

all: build clean

debug: CFLAGS += -g
debug: all

build: $(OBJECTS)
	$(CC) $(LDFLAGS) mincuts $(OBJECTS)

clean:
	rm *.o
