echo "serial"
./lu 8192 128
./lu 8192 64
./lu 8192 32
./lu 8192 16
./lu 8192 8

./lu 8192 120
./lu 8192 100
./lu 8192 80
./lu 8192 60
./lu 8192 40
./lu 8192 20
./lu 8192 10

echo "cpp11"
./lu_cpp 8192 128
./lu_cpp 8192 64
./lu_cpp 8192 32
./lu_cpp 8192 16
./lu_cpp 8192 8

./lu_cpp 8192 120
./lu_cpp 8192 100
./lu_cpp 8192 80
./lu_cpp 8192 60
./lu_cpp 8192 40
./lu_cpp 8192 20
./lu_cpp 8192 10

echo "hpx"
./lu_hpx 8192 128 -t12
./lu_hpx 8192 64 -t12
./lu_hpx 8192 32 -t12
./lu_hpx 8192 16 -t12
./lu_hpx 8192 8 -t12

./lu_hpx 8192 120 -t12
./lu_hpx 8192 100 -t12
./lu_hpx 8192 80 -t12
./lu_hpx 8192 60 -t12
./lu_hpx 8192 40 -t12
./lu_hpx 8192 20 -t12
./lu_hpx 8192 10 -t12

echo "omp"
./lu_omp 8192 128
./lu_omp 8192 64
./lu_omp 8192 32
./lu_omp 8192 16
./lu_omp 8192 8

./lu_omp 8192 120
./lu_omp 8192 100
./lu_omp 8192 80
./lu_omp 8192 60
./lu_omp 8192 40
./lu_omp 8192 20
./lu_omp 8192 10
