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

// NIZKPoK ((g, h); t |h= g^t ) with ROC
// proof is an integer hat_t, and an integer c (lambda_sound bits)
// used in setup - generation of g0

void 
CL_NIZK_RDL_proof(BICYCL::Mpz &hat_t, BICYCL::Mpz &c, BICYCL::QFI h,  BICYCL::QFI g,  BICYCL::Mpz t, BICYCL::Mpz bound_expo, BICYCL::QFIPrecomp g_precomp, BICYCL::CL_HSMqk CL, BICYCL::RandGen &randgen)
{ 
  BICYCL::HashAlgo H{BICYCL::HashAlgo::SHA3_256};
  BICYCL::Mpz bound, tmp_mpz;
  BICYCL::Mpz tilde_t;
  BICYCL::QFI tilde_h;

  BICYCL::QFI tmp_qfi;


  // Commitment

  BICYCL::Mpz::mulby2k (bound, bound_expo, lambda_stat + lambda_sound ); // 2^(2λ_stat + λ_sound) q tilde_s

  //    tilde_h = g^tilde_t , tilde_t < bound
  tilde_t = randgen.random_mpz (bound);
  power_of_g(tilde_h, g, tilde_t, g_precomp, CL);
  
  // Challenge
  
  c = H(h, g, tilde_h);

  // Response

  // hat_t = tilde_t + c t

  hat_t = tilde_t;
  BICYCL::Mpz::addmul(hat_t, c, t);
}


int 
CL_NIZK_RDL_verify(BICYCL::Mpz hat_t, BICYCL::Mpz c, BICYCL::QFI h,  BICYCL::QFI g,  BICYCL::QFIPrecomp g_precomp, BICYCL::CL_HSMqk CL)
{ 
  BICYCL::HashAlgo H{BICYCL::HashAlgo::SHA3_256};
  BICYCL::Mpz verif_c;
  BICYCL::QFI tilde_h;

  BICYCL::QFI tmp_qfi;
  BICYCL::QFI inv_h = h;
  inv_h.neg();  // h^(-1)

  // recompute tilde_h =  g^hat_t * (h^(-1))^c

  power_of_g(tilde_h, g, hat_t, g_precomp, CL);
  CL.Cl_G().nupow(tmp_qfi, inv_h, c);
  CL.Cl_G().nucomp(tilde_h, tilde_h, tmp_qfi);

  // Verif challenge
  
  verif_c = H(h, g, tilde_h);

  return (c == verif_c);
}
