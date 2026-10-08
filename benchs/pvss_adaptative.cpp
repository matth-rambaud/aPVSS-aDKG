/*
 * Efficient Adaptively Secure PVSS from Class Groups and application to DKG 
 * Copyright (C) 2026  
 *                     Guilhem Castagnos 
 *                     Fabien Laguillaumie
 *                     Matthieu Rambaud
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "bicycl.hpp"
#include <thread>
#include <iostream>
#include <fstream>

/* #######################################################################

    
    pvss adaptive

 ####################################################################### */

const size_t lambda = 128; 
const size_t lambda_stat = 40; 
const size_t lambda_sound = 128; 
const size_t q_size = 256;

const unsigned int nb_threads = std::thread::hardware_concurrency();

const unsigned int K = 1; // PVSS parameter -> fixed to 1 now

unsigned int n, t_param;

#include "misc.cpp"
#include "share_decrypt.cpp"

// NIZK proofs

#include "nizk_rdl.cpp"
#include "nizk_rrep.cpp"
#include "nizk_rsh.cpp"
#include "nizk_r2rep.cpp"


/* #######################################################################

main

 ####################################################################### */


int
main (int argc, char* argv[])
{
  if (argc == 3)  {
    n = std::stoi(argv[1]); // Number of P_i
    t_param = std::stoi(argv[2]); // Number of P_i
  }
  else {
    n = 10; // Number of P_i
    t_param = 4; // threshold (t+1 for reconstruction)
  }

  std::ofstream output_file;
  output_file.open ("bench_pvss_"+std::to_string(n)+"_"+std::to_string(t_param)+".csv", std::ios_base::app);

  // print parameters

  std::cout << "Sec Param      : "<< lambda << std::endl;
  std::cout << "Stat Param     : "<< lambda_stat << std::endl;
  std::cout << "Soundness Param: "<< lambda_sound << std::endl;
  std::cout << "Size of q      : "<< q_size << std::endl;
  std::cout << "n              : "<< n << std::endl;
  std::cout << "t              : "<< t_param << std::endl;
  std::cout << std::endl;

  BICYCL::ECGroup EC((BICYCL::SecLevel) 128); // NIST P-256 for lambda= 128
  
  BICYCL::Mpz q = EC.order(); // will be used in setup to create CL group with this q
  // std::cout << q << std::endl;
  
  #include "setup.cpp"

  #include "keygen.cpp"


  std::cout << "************************************ Dist **************************************  " << std::endl;

  std::cout <<std::endl;

  BICYCL::Mpz secret, r_ciphertext;
  std::vector <BICYCL::Mpz> shares (n + 1); // s[i] = P(i) i =0,1,..,n
  std::vector <BICYCL::Mpz> coeffs (t_param + 1); // P = \sum_i=0^t coeff_iX^i
  std::vector <BICYCL::QFI> h_ciphertext (K + 1);
  std::vector <BICYCL::QFI> c_ciphertext (n);

  secret = 42UL;

  start = std::chrono::steady_clock::now();

  r_ciphertext = randgen.random_mpz (bound_expo);

  // generate poly, eval and ciphertext
  share(shares, coeffs, r_ciphertext,  h_ciphertext, c_ciphertext, secret, gg, pk, bound_expo, gg_precomp, pk_precomp, CL, randgen);
  
  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "Share: " << diff.count()*1000 << " ms" <<std::endl;

  BICYCL::Mpz gamma, hat_r;

  start = std::chrono::steady_clock::now();

  CL_NIZK_RSH_proof(gamma, hat_r, pk, gg, h_ciphertext, c_ciphertext, r_ciphertext, bound_expo, gg_precomp, pk_precomp, CL, randgen);

  end = std::chrono::steady_clock::now();
  diff2 = end - start;
  std::cout << "Prove: " << diff2.count()*1000 << " ms" <<std::endl;

  auto dealer_time = diff.count() + diff2.count();
  std::cout << "Dealer Time (share+prove): " << dealer_time * 1000 << " ms" <<std::endl;
 

  
  start = std::chrono::steady_clock::now();

  b = CL_NIZK_RSH_verify(gamma, hat_r, pk, gg, h_ciphertext, c_ciphertext, gg_precomp, pk_precomp, CL);

  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "Verif: " << diff.count()*1000 << " ms" <<std::endl;

  std::cout <<  "Proof ok: " << b << std::endl;

  std::cout <<std::endl;
  std::cout << "************************************ Reconstruction ****************************  " << std::endl;

  std::cout <<std::endl;


  // P0 decryption and proof
  BICYCL::QFI A;
  std::vector <BICYCL::Mpz> s_decrypt (n + 1);

  start = std::chrono::steady_clock::now();

  decrypt(A, s_decrypt[1], h_ciphertext, c_ciphertext[0], x[0], CL);
    
  end = std::chrono::steady_clock::now();
  diff2 = end - start;
  std::cout << "Decrypt: " << diff2.count()*1000 << " ms" <<std::endl;
  
  std::vector <BICYCL::Mpz> c2_rep(n);
  std::vector <std::vector<BICYCL::Mpz>> hat2_x(n, std::vector <BICYCL::Mpz> (K+1));
  
  start = std::chrono::steady_clock::now();
  
  CL_NIZK_2REP_proof(hat2_x[0], c2_rep[0], pk[0], A, gg, h_ciphertext, x[0], bound_expo, gg_precomp, CL, randgen);

  end = std::chrono::steady_clock::now();
  std::chrono::duration<double>  diff3 = end - start;
  std::cout << "Proof: " << diff3.count()*1000 << " ms" <<std::endl;

  auto receiver_time = diff.count() + diff2.count();

  std::cout << "Receiver Time (verif share+ decrypt): " << receiver_time * 1000 << " ms" <<std::endl;


  // P1 ->  Pt
  
  for (unsigned int i = 1; i <= t_param  ; i++) {
    decrypt(A, s_decrypt[i+1], h_ciphertext, c_ciphertext[i], x[i], CL);
    CL_NIZK_2REP_proof(hat2_x[i], c2_rep[i], pk[i], A, gg, h_ciphertext, x[i], bound_expo, gg_precomp, CL, randgen);
  }

  // P0 verif t proof
  
  start = std::chrono::steady_clock::now();

  b = 1;
  
  for (unsigned int i = 1; i <= t_param ; i++) {
    tmp_qfi = CL.power_of_f(s_decrypt[i+1]);
    CL.Cl_G().nucompinv(A, c_ciphertext[i], tmp_qfi);
    
    b &= CL_NIZK_2REP_verify( hat2_x[i], c2_rep[i], pk[i], A, gg, h_ciphertext, gg_precomp,  CL);
  }

  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "Verif: " << diff.count()*1000 << " ms" <<std::endl;

  std::cout << "Proofs ok: " << b << std::endl;


  // Test reconstruction with parts 1 to t+1

  // compute Lagrange's coefficients :
  // for i=1, t+1 : l_i(0) = \prod_{j=1, j≠i}^t+1 -j/(i-j)
  std::vector <BICYCL::Mpz> L (t_param + 2);
  for (long i = 1; i <= t_param + 1 ; i++) {
    L[i] = (1UL);
    for (long j = 1; j <= t_param + 1 ; j++) {
      if (j != i) {
	BICYCL::Mpz::mul(L[i], L[i], j); 
	BICYCL::Mpz::mod_inverse(tmp_mpz, (BICYCL::Mpz) (j - i), CL.q());
	BICYCL::Mpz::mul(L[i], L[i], tmp_mpz);
	BICYCL::Mpz::mod(L[i], L[i], CL.q());
      }
    }
  }

  // s = \sum_i=1^(t+1) s[i] l_i(0)

  tmp_mpz = (0UL);
  for (unsigned int i = 1; i <= t_param + 1 ; i++) {
    BICYCL::Mpz::addmul(tmp_mpz, s_decrypt[i], L[i]);
    BICYCL::Mpz::mod(tmp_mpz, tmp_mpz, CL.q());
  }

  std::cout << "Secret Recontructed: " << tmp_mpz << std::endl;


  output_file << (dealer_time * 1000) <<  ", " << (receiver_time * 1000) << std::endl;
  
  output_file.close();

}
