> :warning: Warning
> This repository is a research project and is **not** production-ready.
> It should **not** be used in real cryptographic projects.

Efficient Adaptively Secure PVSS from Class Groups and application to DKG 
======

This repository contains the implementation in C++ to benchmark the PVSS and DKG with adaptive security from the eprint [2026/xxxx](https://eprint.iacr.org/2026/xxxx).

# Authors

* Guilhem Castagnos 
* Fabien Laguillaumie 
* Matthieu Rambaud 

# Needed Dependencies :

* gmp, openssl

* BICYCL: [GIT](https://gite.lirmm.fr/crypto/bicycl)

# How to build

`mkdir build`

`cd build`

`cmake ..`

`make`

This builds the benchmark in benchs 

# Benchmarks

 `adaptative_pvss n t` : bench PVSS with adaptative security with n parties and threshold t

 `dkg_adaptive n t b` : bench DKG with adaptative security with n parties and threshold t, if b=1 use parallel verification of zk sharing proof, b=0 to use sequential verification 

  `bench.sh` : shell script to launch benchmarks for various parameters     
  
  
# Tuning

* dkg_adaptive uses the 256 bits elliptic curve parametrized in src/bicycl/arith/openssl_wrapper.hpp, line 211. Change to NID_X9_62_prime256v1 for NIST P-256, NID_brainpoolP256t1 for BrainpoolP256t

