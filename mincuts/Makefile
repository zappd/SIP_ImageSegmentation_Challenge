CC = gcc
CFLAGS = -c -Wall
LDFLAGS = -o
DEPENDENCIES = graphCuts.h
OBJECTS = main.o graphCuts.o

%.o: %.c $(DEPENDENCIES)
	$(CC) -c -o $@ $< $(CFLAGS)

all: build clean

debug: CFLAGS += -g
debug: all

build: $(OBJECTS)
	$(CC) $(LDFLAGS) mincuts $(OBJECTS)

clean:
	rm *.o