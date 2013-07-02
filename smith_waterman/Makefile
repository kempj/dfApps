default: hpx-futures

all: omp cpp11 hpx-futures

omp: sw_omp.cpp
	g++ sw_omp.cpp -o sw_omp -fopenmp -O3

cpp11: sw_futures.cpp
	g++ sw_futures.cpp -o sw -std=c++11 -lpthread -O3

hpx-futures: sw_hpx.cpp
	c++ -o sw_hpx sw_hpx.cpp `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=sw_hpx

#hpx-dataflow:
#	c++ -o sw_hpx sw_hpx.cpp `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=sw_hpx
