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

 std::cout << "************************************ Setup *************************************  " << std::endl;

  
  auto start = std::chrono::steady_clock::now();

  std::cout << "******* Class group generation *******" << std::endl;

  // Rand init Bicycl
  BICYCL::Mpz seed;
  BICYCL::RandGen randgen;
  //auto T = std::chrono::system_clock::now();
  seed = (42UL); //static_cast<unsigned long> (T.time_since_epoch().count());
  randgen.set_seed (seed);

BICYCL::CL_HSMqk CL (q, 1, (BICYCL::SecLevel) lambda, randgen, BICYCL::CL_Params{lambda_stat});
  // bitsize q, k=1, seclevel, randgen, cl params for lambda_stat : CL.encrypt_randomness_bound().nbits() will be CL.Cl_DeltaK().class_number_bound().nbits() + lambda_stat -2 
  //  std::cout << CL << std::endl;
  // std::cout << "nb bits of overline s: " << CL.Cl_DeltaK().class_number_bound().nbits() << std::endl; //-> 910-914 for 128 bits secu
  //std::cout << "nb bits of randomness: " << CL.encrypt_randomness_bound().nbits() << std::endl; //-> +lambda_stat -2 bits
  BICYCL::Mpz bound_expo;
  BICYCL::Mpz::mul(bound_expo, CL.encrypt_randomness_bound(), CL.q()); // 2^(λ_stat) q tilde_s 
  BICYCL::Mpz tmp_mpz = randgen.random_mpz(CL.q());
  BICYCL::QFI g,tmp_qfi;
  tmp_qfi = CL.power_of_f(tmp_mpz);
  CL.Cl_G().nucomp(g, CL.h(), tmp_qfi);
  // h returned by setup q-th power -> g := f^x h for DDH-f variant
  // std::cout << g << std::endl;


  BICYCL::QFIPrecomp g_precomp (g, bound_expo.nbits());
  
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double>  diff = end - start;
  auto CLGen = diff.count() ;
  std::cout << "CL.Gen: " << CLGen *1000 << " ms" <<std::endl;

  std::cout << std::endl;


  // bench class group
/*
  BICYCL::QFI myg;

  start = std::chrono::steady_clock::now();
  for (unsigned int i = 0; i < 1000 ; i++) {
    tmp_mpz = randgen.random_mpz(CL.encrypt_randomness_bound());
    CL.power_of_h(myg, tmp_mpz);
  }
  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "Fixed based expo class group: " << diff.count()  << " ms" <<std::endl;




  start = std::chrono::steady_clock::now();
  for (unsigned int i = 0; i < 100 ; i++) {
    tmp_mpz = randgen.random_mpz(CL.encrypt_randomness_bound());
    CL.Cl_G().nupow(tmp_qfi, myg, tmp_mpz, 6);
  }
  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "General expo class group: " << diff.count() *10 << " ms" <<std::endl;

  // bench ec

BICYCL::ECPoint B (EC);
BICYCL::ECPoint tmp_ec (EC);

  start = std::chrono::steady_clock::now();
  for (unsigned int i = 0; i < 100000 ; i++) {
    tmp_mpz = randgen.random_mpz (CL.q());
    B =  BICYCL::ECPoint (EC, BICYCL::BN(tmp_mpz));
  }
  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "Fixed based expo EC: " << diff.count()*1000000 / 100000  << " μs" <<std::endl;

  start = std::chrono::steady_clock::now();
  for (unsigned int i = 0; i < 100000 ; i++) {
    tmp_mpz = randgen.random_mpz (CL.q());
    EC.scal_mul(tmp_ec, BICYCL::BN(tmp_mpz), B);
  }
  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "General expo EC: " << diff.count() *1000000 / 100000  << " μs" <<std::endl;
*/

  std::cout << std::endl;
  std::cout << "******* Interactive generation of g0 *******" << std::endl;

  
  
  auto t = new BICYCL::Mpz[n];
  auto h = new BICYCL::QFI[n];
  auto com = new BICYCL::Mpz[n];
  auto com_r = new BICYCL::Mpz[n];
  std::vector <BICYCL::Mpz> hat_t(n);
  auto c = new BICYCL::Mpz[n];

  BICYCL::HashAlgo H{BICYCL::HashAlgo::SHA3_256}; 

  // P0 computation

  start = std::chrono::steady_clock::now();
    
  t[0] = randgen.random_mpz (bound_expo);
  power_of_g(h[0], g, t[0], g_precomp, CL);  // h0 = g^t0
com_r[0] = randgen.random_mpz ((BICYCL::Mpz) lambda); 
  com[0] = H(h[0], com_r[0]); // commit 

  CL_NIZK_RDL_proof( hat_t[0], c[0], h[0],  g,  t[0], bound_expo, g_precomp, CL, randgen); //prove

  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "Commit and prove: " << diff.count()*1000 << " ms" <<std::endl;

  // P1 -> Pn-1 computation
  for (unsigned int i = 1; i < n ; i++) {
    t[i] = randgen.random_mpz (bound_expo);
    power_of_g(h[i], g, t[i], g_precomp, CL);  // hi = g^ti
    com_r[i] = randgen.random_mpz ((BICYCL::Mpz) lambda); 
    com[i] = H(h[i], com_r[i]); // commit
    CL_NIZK_RDL_proof( hat_t[i], c[i], h[i],  g,  t[i], bound_expo, g_precomp, CL, randgen);
  }

 
  start = std::chrono::steady_clock::now();

  // P0 verification of commit and proof of P1 -> Pn-1

  bool b = 1;

  for (unsigned int i = 1; i < n ; i++) {
    b &= (com[i] == ((BICYCL::Mpz) H(h[i], com_r[i])));
    b &= CL_NIZK_RDL_verify(hat_t[i], c[i], h[i],  g,  g_precomp, CL);
  }

  // g0 = prod h_i
  std::vector <BICYCL::QFI> gg(K+1);

  gg[0] = h[0];
  for (unsigned int i = 1; i < n; i++) {
    CL.Cl_G().nucomp(gg[0], gg[0], h[i]); 
  }

  std::vector <BICYCL::QFIPrecomp> gg_precomp(K+1);
  gg_precomp[0] = BICYCL::QFIPrecomp (gg[0], bound_expo.nbits());


  end = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff2 = end - start;
  std::cout << "Verif: " << diff2.count()*1000 << " ms" <<std::endl;

  std::cout << "Proofs ok: " << b << std::endl;

  auto total_g0 = diff.count() + diff2.count();
  std::cout << "Total: " << total_g0 * 1000 << " ms" <<std::endl;

  std::cout << std::endl;
  std::cout << "******* Interactive generation of g1 *******" << std::endl;


  // P0 computation 
  
// suppose K = 1 here

  BICYCL::Mpz bound;
  BICYCL::Mpz::mulby2k (bound, bound_expo, lambda_sound); // 2^(λ_stat+λ_sound) q tilde_s

  start = std::chrono::steady_clock::now();
  t[0] = randgen.random_mpz (bound);
  power_of_g(h[0], gg[0], t[0], gg_precomp[0], CL);  // h0 = g0^t0
com_r[0] = randgen.random_mpz ((BICYCL::Mpz) lambda); 
  com[0] = H(h[0], com_r[0]); // commit 

  CL_NIZK_RDL_proof( hat_t[0], c[0], h[0],  gg[0],  t[0], bound, gg_precomp[0], CL, randgen); //prove

  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "Commit and prove: " << diff.count()*1000 << " ms" <<std::endl;

  // P1 -> Pn-1 computation
  for (unsigned int i = 1; i < n ; i++) {
    t[i] = randgen.random_mpz (bound);
    power_of_g(h[i], gg[0], t[i], gg_precomp[0], CL);  // hi = g0^ti
    com_r[i] = randgen.random_mpz ((BICYCL::Mpz) lambda); 
    com[i] = H(h[i], com_r[i]); // commit
    CL_NIZK_RDL_proof( hat_t[i], c[i], h[i],  gg[0],  t[i], bound, gg_precomp[0], CL, randgen);
  }

  start = std::chrono::steady_clock::now();

  // P0 verification of proof of P1 -> Pn-1

  b = 1;

for (unsigned int i = 1; i < n ; i++) {
    b &= (com[i] == ((BICYCL::Mpz) H(h[i], com_r[i])));
    b &= CL_NIZK_RDL_verify(hat_t[i], c[i], h[i],  gg[0],  gg_precomp[0], CL);
  }
  
  // g_1 = prod h_i

  gg[1] = h[0];
  for (unsigned int i = 1; i < n; i++) {
    CL.Cl_G().nucomp(gg[1], gg[1], h[i]); 
  }

  gg_precomp[1] = BICYCL::QFIPrecomp (gg[1], bound_expo.nbits());
  


  end = std::chrono::steady_clock::now();
  diff2 = end - start;
  std::cout << "Verif: " << diff2.count()*1000 << " ms" <<std::endl;

  std::cout << "Proofs ok: " << b << std::endl;

  auto total_g1K = diff.count() + diff2.count();

  std::cout << "Total: " << total_g1K * 1000 << " ms" <<std::endl;

  std::cout << std::endl;

  std::cout << "Total Setup: " << (CLGen + total_g0 + total_g1K) * 1000 << " ms" <<std::endl;
  std::cout <<std::endl;
