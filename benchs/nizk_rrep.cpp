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

// NIZKPoK ((pp, pk); ; (η,(x'0,...,x'K)) | pk = \prod_{i=1}^K g_i^{2^-η x'i}) 
// proof is a vector of K+1 integers hat_x, and an integer c (lambda_sound bits)
// used in keygen

void 
CL_NIZK_REP_proof(std::vector <BICYCL::Mpz> &hat_x, BICYCL::Mpz &c, BICYCL::QFI pk, std::vector <BICYCL::QFI> gg,  std::vector <BICYCL::Mpz> x, BICYCL::Mpz bound_expo, std::vector <BICYCL::QFIPrecomp> gg_precomp, BICYCL::CL_HSMqk CL, BICYCL::RandGen &randgen)
{ 
  BICYCL::HashAlgo H{BICYCL::HashAlgo::SHA3_256, lambda_sound/8};
  BICYCL::Mpz bound;
  std::vector <BICYCL::Mpz> tilde_x(K+1);
  BICYCL::QFI tilde_pk, tmp_qfi;


  // Commitment

  BICYCL::Mpz::mulby2k (bound, bound_expo, lambda_stat+lambda_sound); // 2^(2λ_stat+λ_stat) q tilde_s

  for (unsigned int j = 0; j < K+1 ; j++) {
    tilde_x[j] = randgen.random_mpz (bound);
  }

  if (K>2) {
  //  tilde_pk = prod_j=0^K g_j^tilde_x_j with window size 6 (optimal window size for all K)
  CL.Cl_G().nupow(tilde_pk, gg, tilde_x, 6);
  }

  else {
  //  direct expo with precomp and product better for K=1,2
    power_of_g(tilde_pk, gg[0], tilde_x[0], gg_precomp[0], CL);  
    for (unsigned int j = 1; j < K+1 ; j++) {
      power_of_g(tmp_qfi, gg[j], tilde_x[j], gg_precomp[j], CL);  
      CL.Cl_G().nucomp(tilde_pk, tilde_pk, tmp_qfi);
    }
  }

  
  // Challenge
  
  c = H(pk, gg, tilde_pk);

  // Response

  // hat_x[j] = tilde_x[j] + c x[j]
  
  for (unsigned int j = 0; j < K + 1; j++) {
    hat_x[j] = tilde_x[j];
    BICYCL::Mpz::addmul(hat_x[j], x[j], c);
  }
}

int 
CL_NIZK_REP_verify(std::vector <BICYCL::Mpz> hat_x, BICYCL::Mpz c, BICYCL::QFI pk, std::vector <BICYCL::QFI> gg, std::vector <BICYCL::QFIPrecomp> gg_precomp, BICYCL::CL_HSMqk CL)
{ 
  BICYCL::HashAlgo H{BICYCL::HashAlgo::SHA3_256, lambda_sound/8};
  BICYCL::Mpz verif_c;
  BICYCL::QFI tilde_pk, tmp_qfi;

  //recompute  commit tilde_pk = \prod_{j=0}^K g_i^{hat x[j]} pk^{-c}


  if (K>2) {
  //  prod_j=0^K g_j^hat_x_j with window size 6 (optimal window size for all K)
  CL.Cl_G().nupow(tilde_pk, gg, hat_x, 6);
  }

  else {
  //  direct expo with precomp and product better for K=1,2
    power_of_g(tilde_pk, gg[0], hat_x[0], gg_precomp[0], CL);  
    for (unsigned int j = 1; j < K+1 ; j++) {
      power_of_g(tmp_qfi, gg[j], hat_x[j], gg_precomp[j], CL);  
      CL.Cl_G().nucomp(tilde_pk, tilde_pk, tmp_qfi);
    }
  }

  CL.Cl_G().nupow(tmp_qfi, pk, c);
  CL.Cl_G().nucompinv(tilde_pk, tilde_pk, tmp_qfi);


  
  // Verif challenge
  
  verif_c = H(pk, gg, tilde_pk);


  return (c == verif_c) && cg_check(pk, CL);
}
