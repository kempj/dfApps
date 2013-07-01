all:
	c++ -o lu_dist lu-hpx-df-dist.cpp `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=lu_dist

dataflow-local:
	c++ -o lu_df lu-hpx-dataflow.cpp `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=lu_df

futures:
	c++ -o lu_hpx lu-hpx-futures.cpp `pkg-config --cflags --libs hpx_application` -liostreams -DHPX_APPLICATION_NAME=lu_hpx
	
debug:
	g++ -o lu_dist lu-hpx-df-dist.cpp -O0 -rdynamic -fPIC -std=c++0x -Wall -Wno-strict-aliasing -Wsign-promo -Werror=vla -Werror=return-type -fdiagnostics-show-option -Wno-unused-but-set-parameter -Wno-unused-but-set-variable -Werror=uninitialized -Wno-sync-nand -Wno-delete-non-virtual-dtor -pthread -DHPX_HAVE_PARCEL_COALESCING -DBOOST_PARAMETER_MAX_ARITY=7 -DHPX_COROUTINE_ARG_MAX=1 -DHPX_COROUTINE_NO_SEPARATE_CALL_SITES -DHPX_LOG_NO_TSS -DHPX_LOG_NO_TS -DHPX_PLUGIN_SYMBOLS_PREFIX=hpx -DBOOST_BIGINT_HAS_NATIVE_INT64 -DHPX_HAVE_CXX11 -DHPX_HAVE_VERIFY_LOCKS=0 -DHPX_THREAD_GUARD_PAGE=1 -DHPX_THREAD_MAINTAIN_PARENT_REFERENCE=0 -DHPX_THREAD_MAINTAIN_PHASE_INFORMATION=0 -DHPX_THREAD_MAINTAIN_DESCRIPTION=0 -DHPX_THREAD_MAINTAIN_BACKTRACE_ON_SUSPENSION=0 -DHPX_THREAD_MAINTAIN_TARGET_ADDRESS=0 -DHPX_THREAD_MAINTAIN_QUEUE_WAITTIME=1 -DHPX_HAVE_STACKTRACES -DHPX_HAVE_NATIVE_TLS -DNDEBUG -DBOOST_DISABLE_ASSERTS -DHPX_BUILD_TYPE=release -DHPX_UTIL_BIND -DHPX_UTIL_FUNCTION -DHPX_UTIL_TUPLE -DHPX_HAVE_CXX_FLAG_O3 -DHPX_HAVE_C_FLAG_O3 -DHPX_HAVE_CXX_FLAG_RDYNAMIC -DHPX_HAVE_C_FLAG_RDYNAMIC -DHPX_HAVE_CXX_FLAG_FPIC -DHPX_HAVE_C_FLAG_FPIC -DHPX_HAVE_CXX_FLAG_STD_CXX0X -DHPX_HAVE_CXX11_RVALUE_REFERENCES -DHPX_HAVE_CXX11_VARIADIC_TEMPLATES -DHPX_HAVE_CXX11_LAMBDAS -DHPX_HAVE_CXX11_AUTO -DHPX_HAVE_CXX11_DECLTYPE -DBOOST_RESULT_OF_USE_DECLTYPE -DHPX_HAVE_CXX11_STD_UNIQUE_PTR -DBOOST_LOCKFREE_HAVE_CXX11_STD_UNIQUE_PTR -DHPX_HAVE_CXX_FLAG_WALL -DHPX_HAVE_C_FLAG_WALL -DHPX_HAVE_CXX_FLAG_WNO_STRICT_ALIASING -DHPX_HAVE_C_FLAG_WNO_STRICT_ALIASING -DHPX_HAVE_CXX_FLAG_WSIGN_PROMO -DHPX_HAVE_CXX_FLAG_WERROR_VLA -DHPX_HAVE_C_FLAG_WERROR_VLA -DHPX_HAVE_CXX_FLAG_WERROR_RETURN_TYPE -DHPX_HAVE_C_FLAG_WERROR_RETURN_TYPE -DHPX_HAVE_CXX_FLAG_FDIAGNOSTICS_SHOW_OPTION -DHPX_HAVE_C_FLAG_FDIAGNOSTICS_SHOW_OPTION -DHPX_HAVE_CXX_FLAG_WNO_UNUSED_BUT_SET_PARAMETER -DHPX_HAVE_C_FLAG_WNO_UNUSED_BUT_SET_PARAMETER -DHPX_HAVE_CXX_FLAG_WNO_UNUSED_BUT_SET_VARIABLE -DHPX_HAVE_C_FLAG_WNO_UNUSED_BUT_SET_VARIABLE -DHPX_HAVE_CXX_FLAG_WERROR_UNINITIALIZED -DHPX_HAVE_C_FLAG_WERROR_UNINITIALIZED -DHPX_HAVE_CXX_FLAG_WNO_SYNC_NAND -DHPX_HAVE_C_FLAG_WNO_SYNC_NAND -DHPX_HAVE_CXX_FLAG_WNO_DELETE_NON_VIRTUAL_DTOR -DHPX_HAVE_GNU_ALIGNED_16 -DBOOST_ATOMIC_HAVE_GNU_ALIGNED_16 -DHPX_HAVE_GNU_128BIT_INTEGERS -DBOOST_ATOMIC_HAVE_GNU_128BIT_INTEGERS -DHPX_HAVE_RDTSC -D_GNU_SOURCE -DHPX_HAVE_PTHREAD_AFFINITY_NP -DHPX_HAVE_CXX_FLAG_PTHREAD -DHPX_HAVE_C_FLAG_PTHREAD -DHPX_APPLICATION_EXPORTS -DBOOST_ENABLE_ASSERT_HANDLER -I/home/jeremy/hpx/include/hpx/external -I/home/jeremy/boost/include  -Wl,-rpath,:/home/jeremy/hpx/lib/hpx:/home/jeremy/boost/lib -O0 -rdynamic -fPIC -std=c++0x -Wall -Wno-strict-aliasing -Wsign-promo -Werror=vla -Werror=return-type -fdiagnostics-show-option -Wno-unused-but-set-parameter -Wno-unused-but-set-variable -Werror=uninitialized -Wno-sync-nand -Wno-delete-non-virtual-dtor -pthread -L/home/jeremy/hpx/lib/hpx -L/home/jeremy/boost/lib -lhpx -lhpx_init -lhpx_serialization -lboost_date_time -lboost_filesystem -lboost_program_options -lboost_regex -lboost_serialization -lboost_system -lboost_thread -lboost_atomic -lboost_chrono -lhpx -liostreams -DHPX_APPLICATION_NAME=lu_dist

debug-df-local:
	g++ -o lu_df lu-hpx-dataflow.cpp -O0 -g -rdynamic -fPIC -std=c++0x -Wall -Wno-strict-aliasing -Wsign-promo -Werror=vla -Werror=return-type -fdiagnostics-show-option -Wno-unused-but-set-parameter -Wno-unused-but-set-variable -Werror=uninitialized -Wno-sync-nand -Wno-delete-non-virtual-dtor -pthread -DHPX_HAVE_PARCEL_COALESCING -DBOOST_PARAMETER_MAX_ARITY=7 -DHPX_COROUTINE_ARG_MAX=1 -DHPX_COROUTINE_NO_SEPARATE_CALL_SITES -DHPX_LOG_NO_TSS -DHPX_LOG_NO_TS -DHPX_PLUGIN_SYMBOLS_PREFIX=hpx -DBOOST_BIGINT_HAS_NATIVE_INT64 -DHPX_HAVE_CXX11 -DHPX_HAVE_VERIFY_LOCKS=0 -DHPX_THREAD_GUARD_PAGE=1 -DHPX_THREAD_MAINTAIN_PARENT_REFERENCE=0 -DHPX_THREAD_MAINTAIN_PHASE_INFORMATION=0 -DHPX_THREAD_MAINTAIN_DESCRIPTION=0 -DHPX_THREAD_MAINTAIN_BACKTRACE_ON_SUSPENSION=0 -DHPX_THREAD_MAINTAIN_TARGET_ADDRESS=0 -DHPX_THREAD_MAINTAIN_QUEUE_WAITTIME=1 -DHPX_HAVE_STACKTRACES -DHPX_HAVE_NATIVE_TLS -DNDEBUG -DBOOST_DISABLE_ASSERTS -DHPX_BUILD_TYPE=release -DHPX_UTIL_BIND -DHPX_UTIL_FUNCTION -DHPX_UTIL_TUPLE -DHPX_HAVE_CXX_FLAG_O3 -DHPX_HAVE_C_FLAG_O3 -DHPX_HAVE_CXX_FLAG_RDYNAMIC -DHPX_HAVE_C_FLAG_RDYNAMIC -DHPX_HAVE_CXX_FLAG_FPIC -DHPX_HAVE_C_FLAG_FPIC -DHPX_HAVE_CXX_FLAG_STD_CXX0X -DHPX_HAVE_CXX11_RVALUE_REFERENCES -DHPX_HAVE_CXX11_VARIADIC_TEMPLATES -DHPX_HAVE_CXX11_LAMBDAS -DHPX_HAVE_CXX11_AUTO -DHPX_HAVE_CXX11_DECLTYPE -DBOOST_RESULT_OF_USE_DECLTYPE -DHPX_HAVE_CXX11_STD_UNIQUE_PTR -DBOOST_LOCKFREE_HAVE_CXX11_STD_UNIQUE_PTR -DHPX_HAVE_CXX_FLAG_WALL -DHPX_HAVE_C_FLAG_WALL -DHPX_HAVE_CXX_FLAG_WNO_STRICT_ALIASING -DHPX_HAVE_C_FLAG_WNO_STRICT_ALIASING -DHPX_HAVE_CXX_FLAG_WSIGN_PROMO -DHPX_HAVE_CXX_FLAG_WERROR_VLA -DHPX_HAVE_C_FLAG_WERROR_VLA -DHPX_HAVE_CXX_FLAG_WERROR_RETURN_TYPE -DHPX_HAVE_C_FLAG_WERROR_RETURN_TYPE -DHPX_HAVE_CXX_FLAG_FDIAGNOSTICS_SHOW_OPTION -DHPX_HAVE_C_FLAG_FDIAGNOSTICS_SHOW_OPTION -DHPX_HAVE_CXX_FLAG_WNO_UNUSED_BUT_SET_PARAMETER -DHPX_HAVE_C_FLAG_WNO_UNUSED_BUT_SET_PARAMETER -DHPX_HAVE_CXX_FLAG_WNO_UNUSED_BUT_SET_VARIABLE -DHPX_HAVE_C_FLAG_WNO_UNUSED_BUT_SET_VARIABLE -DHPX_HAVE_CXX_FLAG_WERROR_UNINITIALIZED -DHPX_HAVE_C_FLAG_WERROR_UNINITIALIZED -DHPX_HAVE_CXX_FLAG_WNO_SYNC_NAND -DHPX_HAVE_C_FLAG_WNO_SYNC_NAND -DHPX_HAVE_CXX_FLAG_WNO_DELETE_NON_VIRTUAL_DTOR -DHPX_HAVE_GNU_ALIGNED_16 -DBOOST_ATOMIC_HAVE_GNU_ALIGNED_16 -DHPX_HAVE_GNU_128BIT_INTEGERS -DBOOST_ATOMIC_HAVE_GNU_128BIT_INTEGERS -DHPX_HAVE_RDTSC -D_GNU_SOURCE -DHPX_HAVE_PTHREAD_AFFINITY_NP -DHPX_HAVE_CXX_FLAG_PTHREAD -DHPX_HAVE_C_FLAG_PTHREAD -DHPX_APPLICATION_EXPORTS -DBOOST_ENABLE_ASSERT_HANDLER -I/home/jeremy/hpx/include/hpx/external -I/home/jeremy/boost/include  -Wl,-rpath,:/home/jeremy/hpx/lib/hpx:/home/jeremy/boost/lib -O0 -rdynamic -fPIC -std=c++0x -Wall -Wno-strict-aliasing -Wsign-promo -Werror=vla -Werror=return-type -fdiagnostics-show-option -Wno-unused-but-set-parameter -Wno-unused-but-set-variable -Werror=uninitialized -Wno-sync-nand -Wno-delete-non-virtual-dtor -pthread -L/home/jeremy/hpx/lib/hpx -L/home/jeremy/boost/lib -lhpx -lhpx_init -lhpx_serialization -lboost_date_time -lboost_filesystem -lboost_program_options -lboost_regex -lboost_serialization -lboost_system -lboost_thread -lboost_atomic -lboost_chrono -lhpx -liostreams -DHPX_APPLICATION_NAME=lu_df