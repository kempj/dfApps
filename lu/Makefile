
default: dist

all: default hpx-futures omp serial cpp11 dataflow-local

clean: 
	rm -rf lu_cpp lu lu_omp lu_hpx lu_df lu_dist *.o
	
lu_utils_hpx.o: lu_utils.cpp
	g++ -o lu_utils_hpx.o -c lu_utils.cpp -DUSE_HPX `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=lu_utils_hpx.o

lu_utils.o: lu_utils.cpp
	g++ -o lu_utils.o -c lu_utils.cpp -O3 -std=c++11 -lpthread 
	
#g++ -o lu_utils.o -c lu_utils.cpp -std=c++11 -O3 -lpthread 

lu-local.o: lu-local.cpp
	g++ -o lu-local.o -c lu-local.cpp -O3 -std=c++11

dist: lu-hpx-df-dist.cpp
	g++ -o lu_dist lu-hpx-df-dist.cpp `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=lu_dist

dataflow-local: lu-hpx-dataflow.cpp lu_utils_hpx.o lu-local.o
	g++ -o lu_df lu-hpx-dataflow.cpp lu_utils_hpx.o lu-local.o `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=lu_df

df-noglobal: lu-hpx-df-noglobal.cpp lu_utils_hpx.o lu-local.o
	g++ -o lu_df-noglobal lu-hpx-df-noglobal.cpp lu_utils_hpx.o lu-local.o `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=lu_df-noglobal

hpx-futures: lu-hpx-futures.cpp lu_utils_hpx.o lu-local.o
	g++ -o lu_hpx lu-hpx-futures.cpp lu_utils_hpx.o lu-local.o `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=lu_hpx

omp: lu-par-omp-task.c
	g++ lu-par-omp-task.c -o lu_omp -fopenmp -O3

serial:lu-serial.cpp
	g++ lu-serial.cpp -o lu -O3 -fopenmp 

cpp11: lu-futures.cpp lu_utils.o lu-local.o
	g++ lu-futures.cpp lu_utils.o lu-local.o -o lu_cpp -O3 -lpthread -std=c++11
