echo "serial"
./lu 4096 128
./lu 4096 64
./lu 4096 32
./lu 4096 16
./lu 4096 8

./lu 4096 120
./lu 4096 100
./lu 4096 80
./lu 4096 60
./lu 4096 40
./lu 4096 20
./lu 4096 10

echo "cpp11"
./lu_cpp 4096 128
./lu_cpp 4096 64
./lu_cpp 4096 32
./lu_cpp 4096 16
./lu_cpp 4096 8

./lu_cpp 4096 120
./lu_cpp 4096 100
./lu_cpp 4096 80
./lu_cpp 4096 60
./lu_cpp 4096 40
./lu_cpp 4096 20
./lu_cpp 4096 10

echo "hpx"
./lu_hpx 4096 128
./lu_hpx 4096 64
./lu_hpx 4096 32
./lu_hpx 4096 16
./lu_hpx 4096 8

./lu_hpx 4096 120
./lu_hpx 4096 100
./lu_hpx 4096 80
./lu_hpx 4096 60
./lu_hpx 4096 40
./lu_hpx 4096 20
./lu_hpx 4096 10

echo "omp"
./lu_omp 4096 128
./lu_omp 4096 64
./lu_omp 4096 32
./lu_omp 4096 16
./lu_omp 4096 8

./lu_omp 4096 120
./lu_omp 4096 100
./lu_omp 4096 80
./lu_omp 4096 60
./lu_omp 4096 40
./lu_omp 4096 20
./lu_omp 4096 10
