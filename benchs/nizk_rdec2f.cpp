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

// Here we suppose k = 1 

// NIZKPoK ((pp, pk, h, ct, feld);  ((x0,...,xK) ,s) | pk = \prod_{i=1}^K g_i^{xi} ct = f^s \prod_{i=1}^K h_i^{xi}) , feld = s \gprime
// proof is hat_x0, hat_x1, hats, and an integer e (lambda_sound bits)
// used in DKG 2 step

 
void 
CL_NIZK_DEC2f_proof(BICYCL::Mpz &e, BICYCL::Mpz &hat_x0, BICYCL::Mpz &hat_x1, BICYCL::Mpz &hat_s, BICYCL::QFI pk, std::vector <BICYCL::QFI> gg,  std::vector <BICYCL::QFI> h, BICYCL::QFI ct, BICYCL::ECPoint &B, std::vector <BICYCL::Mpz> x, BICYCL::Mpz s, BICYCL::Mpz bound_expo, std::vector <BICYCL::QFIPrecomp> gg_precomp, const BICYCL::ECGroup &EC, BICYCL::CL_HSMqk CL, BICYCL::RandGen &randgen)
{ 
  BICYCL::HashAlgo H{BICYCL::HashAlgo::SHA3_256, lambda_sound/8};
  BICYCL::Mpz bound, tilde_s;
  std::vector <BICYCL::Mpz> tilde_x(K+1);
  BICYCL::QFI tilde_pk, tilde_ct, tmp_qfi;
  BICYCL::ECPoint tilde_B (EC);


  // Commitment

  BICYCL::Mpz::mulby2k (bound, bound_expo, lambda_stat+lambda_sound); // 2^(2λ_stat+λ_stat) q tilde_s

  for (unsigned int j = 0; j < K+1 ; j++) {
    tilde_x[j] = randgen.random_mpz (bound);
  }

  tilde_s = randgen.random_mpz (CL.q());

  //tilde Feld = tilde s P
  
  tilde_B =  BICYCL::ECPoint (EC, BICYCL::BN(tilde_s));
    
  // tilde_pk

  power_of_g(tilde_pk, gg[0], tilde_x[0], gg_precomp[0], CL);  
  power_of_g(tmp_qfi, gg[1], tilde_x[1], gg_precomp[1], CL);  
  CL.Cl_G().nucomp(tilde_pk, tilde_pk, tmp_qfi);

  // tilde_ct

  CL.Cl_G().nupow(tilde_ct, h, tilde_x, 6);
  tmp_qfi = CL.power_of_f(tilde_s);
  CL.Cl_G().nucomp(tilde_ct, tilde_ct, tmp_qfi);
  
    
  // Challenge
  
  e = H(gg, pk, h, ct, ECPointGroupCRefPair(B, EC), tilde_pk, tilde_ct, ECPointGroupCRefPair(tilde_B, EC));

  // Response

  // hat_x[j] = tilde_x[j] + c x[j]
  
  hat_x0 = tilde_x[0];
  BICYCL::Mpz::addmul(hat_x0, x[0], e);
  hat_x1 = tilde_x[1];
  BICYCL::Mpz::addmul(hat_x1, x[1], e);

  // hat_s = tilde_s + s e mod q
  hat_s = tilde_s;
  BICYCL::Mpz::addmul(hat_s, s, e);
  BICYCL::Mpz::mod(hat_s, hat_s, CL.q());
}

int
CL_NIZK_DEC2f_verify(BICYCL::Mpz e, BICYCL::Mpz hat_x0, BICYCL::Mpz hat_x1, BICYCL::Mpz hat_s, BICYCL::QFI pk, std::vector <BICYCL::QFI> gg,  std::vector <BICYCL::QFI> h, BICYCL::QFI ct, BICYCL::ECPoint &B, std::vector <BICYCL::QFIPrecomp> gg_precomp, BICYCL::CL_HSMqk CL)
{ 
  BICYCL::HashAlgo H{BICYCL::HashAlgo::SHA3_256, lambda_sound/8};
  BICYCL::Mpz verif_e;
  BICYCL::QFI tilde_pk, tilde_ct, tmp_qfi;
  

  BICYCL::ECGroup EC((BICYCL::SecLevel) lambda); // can not make work EC in argument with thread, so redeclare here
  BICYCL::ECPoint tilde_B (EC);
    
  bool ret;

  // Check h,ct, pk in hat G

  ret = true;
  for (unsigned int j = 0; j < K + 1; j++) {
    ret &= cg_check(h[j], CL);
  }
  ret &= cg_check(ct, CL);

  // check B in EC

  ret &= EC.is_in_group(B);

  //recompute tilde pk

  power_of_g(tilde_pk, gg[0], hat_x0, gg_precomp[0], CL);  
  power_of_g(tmp_qfi, gg[1], hat_x1, gg_precomp[1], CL);  
  CL.Cl_G().nucomp(tilde_pk, tilde_pk, tmp_qfi);
  CL.Cl_G().nupow(tmp_qfi, pk, e);
  CL.Cl_G().nucompinv(tilde_pk, tilde_pk, tmp_qfi);

  // recompute tilde_ct

  CL.Cl_G().nupow(tilde_ct, h[0], hat_x0);
  CL.Cl_G().nupow(tmp_qfi, h[1], hat_x1);
  CL.Cl_G().nucomp(tilde_ct, tilde_ct, tmp_qfi);
  tmp_qfi = CL.power_of_f(hat_s);
  CL.Cl_G().nucomp(tilde_ct, tilde_ct, tmp_qfi);
  CL.Cl_G().nupow(tmp_qfi, ct, e);
  CL.Cl_G().nucompinv(tilde_ct, tilde_ct, tmp_qfi);

  // recompute tilde_B
  
  BICYCL::ECPoint tmp_ec (EC);
  tilde_B =  BICYCL::ECPoint (EC, BICYCL::BN(hat_s)); //  hat_s P 
  EC.scal_mul(tmp_ec, BICYCL::BN(e), B); //  e.B
  EC.ec_neg(tmp_ec);
  EC.ec_add(tilde_B, tilde_B, tmp_ec); // hat_s P - e.B

  verif_e = H(gg, pk, h, ct, ECPointGroupCRefPair(B, EC), tilde_pk, tilde_ct, ECPointGroupCRefPair(tilde_B, EC));

  ret &= (e == verif_e);
 
  return ret;
}


void
CL_NIZK_DEC2f_verify_multi(std::vector <bool> &res, std::vector <BICYCL::Mpz> e, std::vector <BICYCL::Mpz> hat_x0, std::vector <BICYCL::Mpz> hat_x1, std::vector <BICYCL::Mpz> hat_s, std::vector <BICYCL::QFI> pk, std::vector <BICYCL::QFI> gg, std::vector <BICYCL::QFI> h, std::vector <BICYCL::QFI> c,  std::vector <BICYCL::ECPoint> &B, std::vector <BICYCL::QFIPrecomp> gg_precomp, BICYCL::CL_HSMqk CL, int start, int end) {

  for (int i = start; i < end ; i++) {

    res[i] = CL_NIZK_DEC2f_verify(e[i], hat_x0[i], hat_x1[i], hat_s[i], pk[i], gg,  h, c[i], B[i], gg_precomp, CL);
    
  }
}



