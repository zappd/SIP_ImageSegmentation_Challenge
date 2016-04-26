.PHONY: mincuts docs clean all

all: mincuts docs

mincuts:
	cd mincuts && make

docs:
	cd theory && make pdf

clean:
	cd theory && make clean
	cd mincutes && make clean

