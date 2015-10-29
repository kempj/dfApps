CC=g++

default: hpx-futures

all: omp cpp11 hpx-futures serial

omp: sw_omp.cpp
	$(CC) sw_omp.cpp -o sw_omp -fopenmp -O3

serial: sw_serial.cpp
	$(CC) sw_serial.cpp -o sw -O3

cpp11: sw_futures.cpp
	$(CC) sw_futures.cpp -o sw -std=c++11 -lpthread -O3

hpx-futures: sw_hpx.cpp
	$(CC) -o sw_hpx sw_hpx.cpp `pkg-config --cflags --libs hpx_application` -DHPX_APPLICATION_NAME=sw_hpx

#hpx-dataflow:
#	c++ -o sw_hpx sw_hpx.cpp `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=sw_hpx
