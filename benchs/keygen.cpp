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

  std::cout << "************************************ Keygen ************************************  " << std::endl;

  std::vector <std::vector<BICYCL::Mpz>> x(n, std::vector<BICYCL::Mpz> (K+1)); // nx(K+1) secret key

  std::vector <BICYCL::QFI> pk(n); // n public key
  std::vector <BICYCL::QFIPrecomp> pk_precomp(n);
  
  // P0 computation

  start = std::chrono::steady_clock::now();

  for (unsigned int j = 0; j < K+1 ; j++) {
    x[0][j] = randgen.random_mpz (bound_expo); // secret key
  }

  if (K>2) {
  //  pk = prod_j=0^K g_j^x_j with window size 6 (optimal window size for all K)
  CL.Cl_G().nupow(pk[0], gg, x[0], 6);
  }

  else {
  //  direct expo with precomp and product better for K=1,2
    power_of_g(pk[0], gg[0], x[0][0], gg_precomp[0], CL);  
    for (unsigned int j = 1; j < K+1 ; j++) {
      power_of_g(tmp_qfi, gg[j], x[0][j], gg_precomp[j], CL);  
      CL.Cl_G().nucomp(pk[0], pk[0], tmp_qfi);
    }
  }
  pk_precomp[0] = BICYCL::QFIPrecomp (pk[0], bound_expo.nbits());

  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "pk/sk: " << diff.count()*1000 << " ms" <<std::endl;

  // proof

  start = std::chrono::steady_clock::now();

  std::vector <std::vector<BICYCL::Mpz>> hat_x(n, std::vector <BICYCL::Mpz> (K+1));
  std::vector <BICYCL::Mpz> c_rep(n);

  CL_NIZK_REP_proof(hat_x[0], c_rep[0], pk[0], gg, x[0], bound_expo, gg_precomp, CL, randgen);

  end = std::chrono::steady_clock::now();
  diff2 = end - start;
  std::cout << "Prove: " << diff2.count()*1000 << " ms" <<std::endl;

  std::cout << "Total Keygen : " << (diff.count() + diff2.count()) * 1000 << " ms" <<std::endl;


  // P1 -> Pn-1 computation

  for (unsigned int i = 1; i < n ; i++) {
    for (unsigned int j = 0; j < K+1 ; j++) {
      x[i][j] = randgen.random_mpz (bound_expo);
    }

    if (K>2) {
      //  pk = prod_j=0^K g_j^x_j with window size 6 (optimal window size for all K)
      CL.Cl_G().nupow(pk[i], gg, x[i], 6);
    }

    else {
      //  direct expo with precomp and product better for K=1,2
      power_of_g(pk[i], gg[0], x[i][0], gg_precomp[0], CL);  
      for (unsigned int j = 1; j < K+1 ; j++) {
	power_of_g(tmp_qfi, gg[j], x[i][j], gg_precomp[j], CL);  
	CL.Cl_G().nucomp(pk[i], pk[i], tmp_qfi);
      }
    }
    CL_NIZK_REP_proof(hat_x[i], c_rep[i], pk[i], gg, x[i], bound_expo, gg_precomp, CL, randgen);
    pk_precomp[i] = BICYCL::QFIPrecomp (pk[i], bound_expo.nbits());
  }

  

  
  // time verif. We could do batch verification but with uncompressed proof

  start = std::chrono::steady_clock::now();

  b = 1;
  
  for (unsigned int i = 1; i < n ; i++) {
    b &= CL_NIZK_REP_verify( hat_x[i], c_rep[i], pk[i], gg, gg_precomp, CL);
  }
  
  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "Verif: " << diff.count()*1000 << " ms" <<std::endl;

  std::cout << "Proofs ok: " << b << std::endl;

  std::cout <<std::endl;
