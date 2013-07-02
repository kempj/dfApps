./sgen file_size_50_B 4 2 8192
./sgen file_size_50_A 4 2 8192
./sw 1 file_size_50_A file_size_50_B 1
./sw_omp 1 file_size_50_A file_size_50_B 4096
./sw_hpx 1 file_size_50_A file_size_50_B 4096
