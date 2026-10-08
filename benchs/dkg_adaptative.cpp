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

    
    DKG adaptive

 ####################################################################### */

const size_t lambda = 128; 
const size_t lambda_stat = 40; 
const size_t lambda_sound = 128; 
const size_t q_size = 256;

// window size for W-NAF multi-expo elliptic curve
const int w_ecmulti_exp = 5; // 5 is best on M1-Pro
const  int nb_precomp_w_ecmulti_exp = (1 << (w_ecmulti_exp - 2));


const unsigned int nb_threads = std::thread::hardware_concurrency();

const unsigned int K = 1; // PVSS parameter -> fixed to 1 now

unsigned int n, t_param;

bool parallel_verif;

#include "misc.cpp"
#include "share_decrypt.cpp"

// NIZK proofs

#include "nizk_rdl.cpp"
#include "nizk_rrep.cpp"
#include "ec_multi_exp.cpp"
#include "nizk_rfsh.cpp"
#include "nizk_rdec2f.cpp"




/* #######################################################################

main

 ####################################################################### */


int
main (int argc, char* argv[])
{
  if (argc == 4)  {
    n = std::stoi(argv[1]); // Number of P_i
    t_param = std::stoi(argv[2]); // threshold (t+1 for reconstruction)
    parallel_verif = std::stoi(argv[3]); // 1 if parallel verif, 0 otherwise
  }
  else {
    n = 10; // Number of P_i
    t_param = 4; // threshold (t+1 for reconstruction)
    parallel_verif = 1;
  }

  std::ofstream output_file;
  output_file.open ("bench_dkg_"+std::to_string(n)+"_"+std::to_string(t_param)+"_"+std::to_string(parallel_verif)+".csv", std::ios_base::app);

  // print parameters
  
  std::cout << "Sec Param      : "<< lambda << std::endl;
  std::cout << "Stat Param     : "<< lambda_stat << std::endl;
  std::cout << "Soundness Param: "<< lambda_sound << std::endl;
  std::cout << "Size of q      : "<< q_size << std::endl;
  std::cout << "n              : "<< n << std::endl;
  std::cout << "t              : "<< t_param << std::endl;
  std::cout << std::endl;

  const BICYCL::ECGroup EC((BICYCL::SecLevel) lambda); // NIST P-256 for lambda= 128

  BICYCL::Mpz q = EC.order(); // will be used in setup to create CL group with this q

  std::cout << "q:"<< q << std::endl;

  #include "setup.cpp"

  #include "keygen.cpp"

  std::cout << "************************************ Dist **************************************  " << std::endl;

  std::cout <<std::endl;

  // ################### Step 1 Fdist
  
  BICYCL::Mpz secret; // s_i player i no need to keep that
  std::vector <std::vector<BICYCL::Mpz>> shares (n, std::vector<BICYCL::Mpz> (n + 1)); // shares[i][j] = Qi(j) i =0,1,..,n-1 j=1,..n
  std::vector <std::vector<BICYCL::Mpz>> coeffs (n, std::vector<BICYCL::Mpz> (t_param + 1)); // Q_i = \sum_k=0^t coeff[i][k] X^k  i=0,1,...n-1
  std::vector<BICYCL::Mpz> r_ciphertext (n); 
  std::vector <std::vector <BICYCL::QFI>> h_ciphertext(n,  std::vector <BICYCL::QFI> (K + 1));
  std::vector <std::vector <BICYCL::QFI>> c_ciphertext(n,  std::vector <BICYCL::QFI> (n));
  std::vector <BICYCL::ECPoint> A_commit;

  std::vector <BICYCL::Mpz> nizkfsh_gamma (n);
  std::vector <BICYCL::Mpz> nizkfsh_hat_r (n);
  std::vector <BICYCL::Mpz> nizkfsh_hat_s (n);


  
  // P0 computation

  start = std::chrono::steady_clock::now();

  secret = randgen.random_mpz(CL.q());

  r_ciphertext[0] = randgen.random_mpz (bound_expo);

  share(shares[0], coeffs[0], r_ciphertext[0],  h_ciphertext[0], c_ciphertext[0], secret, gg, pk, bound_expo, gg_precomp, pk_precomp, CL, randgen);

  A_commit.push_back(BICYCL::ECPoint (EC, BICYCL::BN(secret)));

  CL_NIZK_RfSH_proof(nizkfsh_gamma[0], nizkfsh_hat_r[0], nizkfsh_hat_s[0], pk, gg, h_ciphertext[0], c_ciphertext[0], A_commit[0], secret, r_ciphertext[0], bound_expo, gg_precomp, pk_precomp, EC, CL, randgen);
  
  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "fDist (Share + Proof fSh): " << diff.count()*1000 << " ms" <<std::endl;

  // P1 -> Pn computation

  for (unsigned int i = 1; i < n ; i++) {

    secret = randgen.random_mpz(CL.q());
 
    r_ciphertext[i] = randgen.random_mpz (bound_expo);
  
    share(shares[i], coeffs[i], r_ciphertext[i],  h_ciphertext[i], c_ciphertext[i], secret, gg, pk, bound_expo, gg_precomp, pk_precomp, CL, randgen);

    A_commit.push_back(BICYCL::ECPoint (EC, BICYCL::BN(coeffs[i][0])));
       
     
    CL_NIZK_RfSH_proof(nizkfsh_gamma[i], nizkfsh_hat_r[i], nizkfsh_hat_s[i], pk, gg, h_ciphertext[i], c_ciphertext[i], A_commit[i], secret, r_ciphertext[i], bound_expo, gg_precomp, pk_precomp, EC, CL, randgen);
     

  }

    // ################### Step 2
  
  // verification of t proofs by P0 VerifySharing

  
  start = std::chrono::steady_clock::now();

  bool ok = 1;

  if (!parallel_verif) {
    
    // without threads
    
    for (unsigned int i = 1; i < t_param + 1; i++) {
       ok &= CL_NIZK_RfSH_verify(nizkfsh_gamma[i], nizkfsh_hat_r[i], nizkfsh_hat_s[i], pk, gg, h_ciphertext[i], c_ciphertext[i], A_commit[i], gg_precomp, pk_precomp, CL);
    }
  }

  else {
  
    // with threads

    std::vector <bool> res (n);
    std::cout << "Parallel ";
    auto th = new std::thread[nb_threads];
    unsigned int nb = t_param; 
    unsigned int local_nb_threads;
    unsigned int quotient, remainder, first, last;

    if (nb < nb_threads) {
      local_nb_threads = nb;
      quotient = 1;
      remainder = 0;
    }
    else {
      local_nb_threads = nb_threads;
      quotient = nb/nb_threads;
      remainder = nb % nb_threads;
    }
    
    first = 1;
    for (unsigned int tt = 0; tt < local_nb_threads - 1; tt++) {
      last = first + quotient + (tt < remainder ? 1 : 0);
      th[tt] = std::thread(CL_NIZK_RfSH_verify_multi, std::ref(res), nizkfsh_gamma, nizkfsh_hat_r, nizkfsh_hat_s, pk, gg, h_ciphertext, c_ciphertext, std::ref(A_commit), gg_precomp, pk_precomp, CL, first, last);
      first = last;
    }
    
    CL_NIZK_RfSH_verify_multi(std::ref(res), nizkfsh_gamma, nizkfsh_hat_r, nizkfsh_hat_s, pk, gg, h_ciphertext, c_ciphertext, A_commit, gg_precomp, pk_precomp, CL, first, nb + 1);

    for (unsigned int tt = 0; tt < local_nb_threads-1; tt++) {
      th[tt].join();
    }

    for (unsigned int i = 1; i < t_param + 1; i++) {
      ok &= res[i];
    }
  }
    
  end = std::chrono::steady_clock::now();
  diff2 = end - start;
  std::cout << "Verif fSh: " << diff2.count()*1000 << " ms " << std::endl;
  std::cout << "Proofs ok: " << ok << std::endl;


   
  //  for U = {t+1 first players}

  //tpk = sum_U Feld_j

  start = std::chrono::steady_clock::now();
    
  BICYCL::ECPoint tpk (EC,  A_commit[0]);
  for (unsigned int i = 1; i < t_param + 1  ; i++) {
    EC.ec_add(tpk, tpk, A_commit[i]);
  }

  
  std::vector <BICYCL::QFI> h_tsk (K + 1); 
  std::vector <BICYCL::QFI> c_tsk (n);


  // P0: aggregate ciphertexts for all players in U
  // (\oplus_{j \in U} h_j,0, \oplus h_j,1) (\oplus ct_j,0,..., \oplus ct_j,n)
  
  for (unsigned int j = 0; j < K + 1; j++) {
    h_tsk[j] =  h_ciphertext[0][j];
  }
  for (unsigned int j = 0; j < n; j++) {
    c_tsk[j] =  c_ciphertext[0][j];
  }
  
  for (unsigned int i = 1; i < t_param + 1; i++) {
    for (unsigned int j = 0; j < K + 1; j++) {
      CL.Cl_G().nucomp(h_tsk[j], h_tsk[j], h_ciphertext[i][j]);
    }
    for (unsigned int j = 0; j < n; j++) {
      CL.Cl_G().nucomp(c_tsk[j], c_tsk[j], c_ciphertext[i][j]);
    }
  }

  std::vector <BICYCL::Mpz> s (t_param + 1);
  std::vector <BICYCL::ECPoint> B_commit;
    
  // P0: Dec2Feldman : decrypt its ciphertext (htsk, c_tsk_0), commit and prove
  // Here we suppose K = 1


  decrypt(tmp_qfi, s[0], h_tsk, c_tsk[0], x[0], CL);
  B_commit.push_back(BICYCL::ECPoint (EC, BICYCL::BN(s[0])));

  std::vector <BICYCL::Mpz> nizkdec2f_e (t_param + 1);
  std::vector <BICYCL::Mpz> nizkdec2f_hat_x0 (t_param + 1);
  std::vector <BICYCL::Mpz> nizkdec2f_hat_x1 (t_param + 1);
  std::vector <BICYCL::Mpz> nizkdec2f_hat_s (t_param + 1);
  
  CL_NIZK_DEC2f_proof(nizkdec2f_e[0], nizkdec2f_hat_x0[0], nizkdec2f_hat_x1[0], nizkdec2f_hat_s[0], pk[0], gg, h_tsk, c_tsk[0], B_commit[0], x[0], s[0], bound_expo, gg_precomp, EC, CL, randgen);

  end = std::chrono::steady_clock::now();
  std::chrono::duration<double>   diff3 = end - start;
  std::cout << "Dec2Feldman: " << diff3.count()*1000 << " ms " << std::endl;



  // P1 -> Pn computation

  for (unsigned int i = 1; i < t_param + 1 ; i++) {
    decrypt(tmp_qfi, s[i], h_tsk, c_tsk[i], x[i], CL);
    B_commit.push_back(BICYCL::ECPoint (EC, BICYCL::BN(s[i])));
    CL_NIZK_DEC2f_proof(nizkdec2f_e[i], nizkdec2f_hat_x0[i], nizkdec2f_hat_x1[i], nizkdec2f_hat_s[i], pk[i], gg, h_tsk, c_tsk[i], B_commit[i], x[i], s[i], bound_expo, gg_precomp, EC, CL, randgen);
  }

  // P0 verif t proof

  start = std::chrono::steady_clock::now();

  ok = 1;

  if (!parallel_verif) {
    
    // without threads
    
    for (unsigned int i = 1; i < t_param + 1; i++) {
      ok &= CL_NIZK_DEC2f_verify(nizkdec2f_e[i], nizkdec2f_hat_x0[i], nizkdec2f_hat_x1[i], nizkdec2f_hat_s[i], pk[i], gg,  h_tsk, c_tsk[i], B_commit[i], gg_precomp, CL);
    }
  }

  else {
  
    // with threads
    std::vector <bool> res2 (n);
    std::cout << "Parallel ";
    auto th = new std::thread[nb_threads];
    unsigned int nb = t_param; 
    unsigned int local_nb_threads;
    unsigned int quotient, remainder, first, last;

    if (nb < nb_threads) {
      local_nb_threads = nb;
      quotient = 1;
      remainder = 0;
    }
    else {
      local_nb_threads = nb_threads;
      quotient = nb/nb_threads;
      remainder = nb % nb_threads;
    }
    
    first = 1;
    for (unsigned int tt = 0; tt < local_nb_threads - 1; tt++) {
      last = first + quotient + (tt < remainder ? 1 : 0);
      th[tt] = std::thread(CL_NIZK_DEC2f_verify_multi, std::ref(res2), nizkdec2f_e, nizkdec2f_hat_x0, nizkdec2f_hat_x1, nizkdec2f_hat_s, pk, gg,  h_tsk, c_tsk, std::ref(B_commit), gg_precomp, CL,first,last);
      first = last;
    }

    CL_NIZK_DEC2f_verify_multi(std::ref(res2), nizkdec2f_e, nizkdec2f_hat_x0, nizkdec2f_hat_x1, nizkdec2f_hat_s, pk, gg,  h_tsk, c_tsk, B_commit, gg_precomp, CL,first,nb+1);
  
    for (unsigned int tt = 0; tt < local_nb_threads-1; tt++) {
      th[tt].join();
    }

    for (unsigned int i = 1; i < t_param ; i++) {
      ok &= res2[i];
    }
  }
    
  end = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff4 = end - start;
  std::cout << "Verif Dec2f: " << diff4.count()*1000 << " ms " << std::endl;
  std::cout << "Proofs ok: " << ok << std::endl;


  auto total_time = diff.count() + diff2.count() + diff3.count() + diff4.count();

  std::cout << "Total Time (fDist + Verify fSh + Dec2Felman + Verif Dec2f " << total_time * 1000 << " ms" <<std::endl;

  // Test reconstruction with parts 1 to t+1

  // compute Lagrange's coefficients :
  // for i=1, t+1 : l_(i-1)(0) = \prod_{j=1, j≠i}^t+1 -j/(i-j)
  std::vector <BICYCL::Mpz> L (t_param + 1);
  for (long i = 1; i <= t_param + 1 ; i++) {
    L[i-1] = (1UL);
    for (long j = 1; j <= t_param + 1 ; j++) {
      if (j != i) {
	BICYCL::Mpz::mul(L[i-1], L[i-1], j); 
	BICYCL::Mpz::mod_inverse(tmp_mpz, (BICYCL::Mpz) (j - i), CL.q());
	BICYCL::Mpz::mul(L[i-1], L[i-1], tmp_mpz);
	BICYCL::Mpz::mod(L[i-1], L[i-1], CL.q());
      }
    }
  }

  // s = \sum_i=0^(t) s[i] l_i(0)

  tmp_mpz = (0UL);
  for (unsigned int i = 0; i < t_param + 1 ; i++) {
    BICYCL::Mpz::addmul(tmp_mpz, s[i], L[i]);
    BICYCL::Mpz::mod(tmp_mpz, tmp_mpz, CL.q());
  }

  // test tpk

  BICYCL::ECPoint test_tpk (EC,  BICYCL::BN(tmp_mpz));
  EC.ec_neg(test_tpk);
  EC.ec_add(test_tpk, tpk ,test_tpk);
  std::cout << "Private shares ok: " << (EC.is_at_infinity(test_tpk)) << std::endl;

  // test tsk
  EC_multi_exp(test_tpk, B_commit,  L, t_param + 1, EC);
  EC.ec_neg(test_tpk);
  EC.ec_add(test_tpk, tpk ,test_tpk);
  std::cout << "Public shares ok: " << (EC.is_at_infinity(test_tpk)) << std::endl;
  
  output_file << (total_time * 1000) << std::endl;
  
  output_file.close();
}
