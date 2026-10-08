include_directories (${BICYCL_INCLUDE_DIRS} ${GMP_INCLUDE_DIRS})
link_libraries (${GMP_LIBRARIES} ${GMP_CPP_LIBRARY} ${OPENSSL_CRYPTO_LIBRARY})


### List of benchmarks
set (BENCHS_LIST pvss_adaptative dkg_adaptative)

### Automatically create target and test from TESTS_LIST
foreach (bench ${BENCHS_LIST})
  add_executable (${bench} ${bench}.cpp)

endforeach ()
