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

// NIZKPoK (\pp, (\pk_1,\dots,\pk_n), (h_0,h_1, \dots, h_K,  c_1,\dots,c_n), \Feld  );
//    ( P(X), r ) \; | \;
//   P(X)=\sum a_iX^i \in \F_q[X]_{\leqslant t},  r \in \Z, (h_0,h_1, \dots, h_K, c_1,\dots,c_n) = (g_0^r, g_1^r,\dots,g_K^r,     f^{P(\alpha_1)}\, \pk_1^r, \dots , f^{P(\alpha_n)}\, \pk_n^r), P(0) \gprime
// proof is  * gamma challenge of lambda_sound bits
//           * hat r : integer of < 2^2lambda_stat+lmanda_sound q \tilde s
//           * hat s : integer < q
// used in sharing DKG

using ECPointGroupCRefPair = std::tuple<const BICYCL::ECPoint &, const BICYCL::ECGroup &>;


void 
CL_NIZK_RfSH_proof(BICYCL::Mpz &gamma, BICYCL::Mpz &hat_r, BICYCL::Mpz &hat_s, std::vector <BICYCL::QFI> pk, std::vector <BICYCL::QFI> gg, std::vector <BICYCL::QFI> h, std::vector <BICYCL::QFI> c, BICYCL::ECPoint &A_commit, BICYCL::Mpz secret, BICYCL::Mpz r, BICYCL::Mpz bound_expo, std::vector <BICYCL::QFIPrecomp> gg_precomp, std::vector <BICYCL::QFIPrecomp> pk_precomp, const BICYCL::ECGroup &EC, BICYCL::CL_HSMqk CL, BICYCL::RandGen &randgen)
{
  BICYCL::HashAlgo H1{BICYCL::HashAlgo::SHAKE128, std::max(32UL, q_size/8* (n-t_param-1) + lambda_sound/8 * n)};
  BICYCL::HashAlgo H2{BICYCL::HashAlgo::SHA3_256, lambda_sound/8};
  BICYCL::Mpz bound, tilde_r, tilde_s, tmp_mpz;
  std::vector <BICYCL::Mpz> m_star(n-t_param-1);
  std::vector <BICYCL::Mpz> e(n);
  std::vector <BICYCL::Mpz> w(n+1);
  std::vector <BICYCL::Mpz> u(n);
    
  std::vector <BICYCL::QFI> tilde_h(K+1);
  BICYCL::QFI U, tilde_V, tmp_qfi;
  BICYCL::ECPoint tilde_A_commit (EC);

  // First challenge w'

  std::vector <unsigned char> tmp_hash;
  std::vector <unsigned char> tmp_hash_extract;
  std::vector<unsigned char>::const_iterator first, last;
  
  tmp_hash = H1(gg, pk, h, c, ECPointGroupCRefPair(A_commit, EC)); // give (q_size/8* (n-t_param-1) + lambda_sound/8 * n)  char
  // the first q_size/8* (n-t_param-1) are used to build the polynomial m_star
  // the following  lambda_sound/8 * n to generate n challenge e

  for (unsigned int j = 0; j < n - t_param - 1; j++) {
    first = tmp_hash.begin() + (q_size/8) * j;
    last = tmp_hash.begin() + (q_size/8) * (j+1);
    tmp_hash_extract = std::vector <unsigned char> (first,last);
    m_star[j] = (BICYCL::Mpz) tmp_hash_extract;
  }

  for (unsigned int j = 0; j < n; j++) {
    first = tmp_hash.begin() + (q_size/8)*(n - t_param -1) + (lambda_sound/8) * j;
    last = tmp_hash.begin() + (q_size/8)*(n - t_param -1)  + (lambda_sound/8) * (j+1);
    tmp_hash_extract = std::vector <unsigned char> (first,last);
    e[j] = (BICYCL::Mpz) tmp_hash_extract;
  }

  // compute v_i coefficients (from 1 to n to avoid headache)
  // for i=1, n : vi = \prod_{j=1, j≠i}^n 1/(i-j)
  std::vector <BICYCL::Mpz> v (n + 1);
  for (long i = 1; i <= n ; i++) {

    v[i] = (1UL);
    for (long j = 1; j <= n ; j++) {
      if (j != i) {
	BICYCL::Mpz::mod_inverse(tmp_mpz, (BICYCL::Mpz) (i - j), CL.q());
	BICYCL::Mpz::mul(v[i], v[i], tmp_mpz);
	BICYCL::Mpz::mod(v[i], v[i], CL.q());
      }
    }
  }

  // compute l_i coefficients (from 1 to n to avoid headache)
  // for i=1, t+1 : vi = \prod_{j=1, j≠i}^(t+1) j/(j-i)
  std::vector <BICYCL::Mpz> L (n + 1);

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

  for (long i = t_param + 2; i <= n ; i++) {
    L[i] = (0UL);
  }

  
  // compute w_i (still from 1 to n)

  // Compute evalutations with Horner's method in i=1,...,n of M_star
    
  for (unsigned int i = 1; i < n + 1; i++) {
    w[i] = m_star[n - t_param - 2];
    for (int k = n - t_param - 3; k >= 0; k--)
      {
	BICYCL::Mpz::mul(w[i], w[i], i);
	BICYCL::Mpz::add(w[i], w[i], m_star[k]);
	BICYCL::Mpz::mod(w[i], w[i], CL.q());
      }

    BICYCL::Mpz::mul(w[i], w[i], v[i]); // multiply eval by v[i] mod q
    BICYCL::Mpz::mod(w[i], w[i], CL.q());
  }

  // now compute u (from 0 to n-1) in Z
  for (unsigned int i = 0; i < n ; i++) { 
    BICYCL::Mpz::mul(tmp_mpz, e[i], CL.q());
    BICYCL::Mpz::add(u[i], w[i+1], tmp_mpz);
    BICYCL::Mpz::add(u[i], L[i+1], u[i]);
  }

  // Compute U

  if (n> 3) {
    CL.Cl_G().nupow(U, pk, u, 6); // multi expo, size 6
  }

  else {
  //  direct expo with precomp and product better for n=2,3
    power_of_g(U, pk[0], u[0], pk_precomp[0], CL);  
    for (unsigned int j = 1; j < n ; j++) {
      power_of_g(tmp_qfi, pk[j], u[j], pk_precomp[j], CL);  
      CL.Cl_G().nucomp(U, U, tmp_qfi);
    }
  }

  
  // Compression done, sigma protocol for CLDL

  // Commitment

  BICYCL::Mpz::mulby2k (bound, bound_expo, lambda_stat+lambda_sound); // 2^(2λ_stat+λ_sound) q tilde_s

  tilde_r = randgen.random_mpz (bound);
  tilde_s = randgen.random_mpz (CL.q());

  // Compute tilde h, tilde V, tilde Feld

  // without threads 
  // g0^r,...,gK^r
  
  for (unsigned int j = 0; j < K + 1; j++) {
    power_of_g(tilde_h[j], gg[j], tilde_r, gg_precomp[j], CL);  
  }

  // tilde V = U^tilde r, f^tilde s
  CL.Cl_G().nupow(tilde_V, U, tilde_r);
  tmp_qfi = CL.power_of_f(tilde_s);
  CL.Cl_G().nucomp(tilde_V, tilde_V, tmp_qfi);

  //tilde Feld = tilde s P

  tilde_A_commit =  BICYCL::ECPoint (EC, BICYCL::BN(tilde_s));
  
  // Challenge
  gamma = H2(gg, pk, h, c, ECPointGroupCRefPair(A_commit, EC), tilde_h, tilde_V, ECPointGroupCRefPair(tilde_A_commit, EC)); 

  // Response

  // hat_r = tilde_r + r gamma
  // hat_s = tilde_s + s gamma mod q

  hat_r = tilde_r;
  BICYCL::Mpz::addmul(hat_r, r, gamma);

  hat_s = tilde_s;
  BICYCL::Mpz::addmul(hat_s, secret, gamma);
  BICYCL::Mpz::mod(hat_s, hat_s, CL.q());
}

int CL_NIZK_RfSH_verify(BICYCL::Mpz gamma, BICYCL::Mpz hat_r, BICYCL::Mpz hat_s, std::vector <BICYCL::QFI> pk, std::vector <BICYCL::QFI> gg, std::vector <BICYCL::QFI> h, std::vector <BICYCL::QFI> c, BICYCL::ECPoint &A_commit, std::vector <BICYCL::QFIPrecomp> gg_precomp, std::vector <BICYCL::QFIPrecomp> pk_precomp, BICYCL::CL_HSMqk CL)
{

  BICYCL::ECGroup EC((BICYCL::SecLevel) lambda); // can not make work EC in argument with thread, so redeclare here

  bool ret;

  // Check h,c in hat G

  ret = true;
  for (unsigned int j = 0; j < K + 1; j++) {
    ret &= cg_check(h[j], CL);
  }
  for (unsigned int j = 0; j < n; j++) {
    ret &= cg_check(c[j], CL);
  }

  // check A in EC

  ret &= EC.is_in_group(A_commit);
  
  BICYCL::HashAlgo H1{BICYCL::HashAlgo::SHAKE128, std::max(32UL, q_size/8* (n-t_param-1) + lambda_sound/8 * n)};
  BICYCL::HashAlgo H2{BICYCL::HashAlgo::SHA3_256, lambda_sound/8};
  BICYCL::Mpz bound, tilde_r, tilde_s, tmp_mpz;
  std::vector <BICYCL::Mpz> m_star(n-t_param-1);
  std::vector <BICYCL::Mpz> e(n);
  std::vector <BICYCL::Mpz> w(n+1);
  std::vector <BICYCL::Mpz> u(n);
    
  std::vector <BICYCL::QFI> tilde_h(K+1);
  BICYCL::QFI U, V, tilde_V, tmp_qfi;
  BICYCL::ECPoint tilde_A_commit (EC);

  // First challenge w'

  std::vector <unsigned char> tmp_hash;
  std::vector <unsigned char> tmp_hash_extract;
  std::vector<unsigned char>::const_iterator first, last;
  
  tmp_hash = H1(gg, pk, h, c, ECPointGroupCRefPair(A_commit, EC)); // give (q_size/8* (n-t_param-1) + lambda_sound/8 * n)  char
  // the first q_size/8* (n-t_param-1) are used to build the polynomial m_star
  // the following  lambda_sound/8 * n to generate n challenge e

  for (unsigned int j = 0; j < n - t_param - 1; j++) {
    first = tmp_hash.begin() + (q_size/8) * j;
    last = tmp_hash.begin() + (q_size/8) * (j+1);
    tmp_hash_extract = std::vector <unsigned char> (first,last);
    m_star[j] = (BICYCL::Mpz) tmp_hash_extract;
  }

  for (unsigned int j = 0; j < n; j++) {
    first = tmp_hash.begin() + (q_size/8)*(n - t_param -1) + (lambda_sound/8) * j;
    last = tmp_hash.begin() + (q_size/8)*(n - t_param -1)  + (lambda_sound/8) * (j+1);
    tmp_hash_extract = std::vector <unsigned char> (first,last);
    e[j] = (BICYCL::Mpz) tmp_hash_extract;
  }

  // compute v_i coefficients (from 1 to n to avoid headache)
  // for i=1, n : vi = \prod_{j=1, j≠i}^n 1/(i-j)
  std::vector <BICYCL::Mpz> v (n + 1);
  for (long i = 1; i <= n ; i++) {

    v[i] = (1UL);
    for (long j = 1; j <= n ; j++) {
      if (j != i) {
	BICYCL::Mpz::mod_inverse(tmp_mpz, (BICYCL::Mpz) (i - j), CL.q());
	BICYCL::Mpz::mul(v[i], v[i], tmp_mpz);
	BICYCL::Mpz::mod(v[i], v[i], CL.q());
      }
    }
  }

    // compute l_i coefficients (from 1 to n to avoid headache)
  // for i=1, t+1 : vi = \prod_{j=1, j≠i}^(t+1) j/(j-i)

    std::vector <BICYCL::Mpz> L (n + 1);

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

  for (long i = t_param + 2; i <= n ; i++) {
    L[i] = (0UL);
  }

  // compute w_i (still from 1 to n)

  // Compute evalutations with Horner's method in i=1,...,n of M_star
    
  for (unsigned int i = 1; i < n + 1; i++) {
    w[i] = m_star[n - t_param - 2];
    for (int k = n - t_param - 3; k >= 0; k--)
      {
	BICYCL::Mpz::mul(w[i], w[i], i);
	BICYCL::Mpz::add(w[i], w[i], m_star[k]);
	BICYCL::Mpz::mod(w[i], w[i], CL.q());
      }

    BICYCL::Mpz::mul(w[i], w[i], v[i]); // multiply eval by v[i] mod q
    BICYCL::Mpz::mod(w[i], w[i], CL.q());
  }

  // now compute u (from 0 to n-1) in Z
  for (unsigned int i = 0; i < n ; i++) { 
    BICYCL::Mpz::mul(tmp_mpz, e[i], CL.q());
    BICYCL::Mpz::add(u[i], w[i+1], tmp_mpz);
    BICYCL::Mpz::add(u[i], L[i+1], u[i]);
  }

  // Compute U

  if (n> 3) {
    CL.Cl_G().nupow(U, pk, u, 6); // multi expo, size 6
  }

  else {
  //  direct expo with precomp and product better for n=2,3
    power_of_g(U, pk[0], u[0], pk_precomp[0], CL);  
    for (unsigned int j = 1; j < n ; j++) {
      power_of_g(tmp_qfi, pk[j], u[j], pk_precomp[j], CL);  
      CL.Cl_G().nucomp(U, U, tmp_qfi);
    }
  }

  // Compute V

  CL.Cl_G().nupow(V, c, u, 6); // multi expo, size 6. (no precomp for c)


  // Compression done, sigma protocol for CLDL

  //recompute tilde_h_j = g_j^{hat r} h_j^{-gamma}

  // without threads

  for (unsigned int j = 0; j < K + 1; j++) {
    power_of_g(tilde_h[j], gg[j], hat_r, gg_precomp[j], CL);
    CL.Cl_G().nupow(tmp_qfi, h[j], gamma);
    CL.Cl_G().nucompinv(tilde_h[j], tilde_h[j], tmp_qfi);
  }

  // recompute tilde V = U^hat_r f^hat s V^{-gamma}

  
  CL.Cl_G().nupow(tilde_V, U, hat_r);
  tmp_qfi = CL.power_of_f(hat_s);
  CL.Cl_G().nucomp(tilde_V, tilde_V, tmp_qfi);
  CL.Cl_G().nupow(tmp_qfi, V, gamma);
  CL.Cl_G().nucompinv(tilde_V, tilde_V, tmp_qfi);

  // recompute tilde Feld = hat_s P - gamma B

  BICYCL::ECPoint tmp_ec (EC);
  EC.scal_mul(tmp_ec, BICYCL::BN(gamma), A_commit); // gamma Feld
  tilde_A_commit =  BICYCL::ECPoint (EC, BICYCL::BN(hat_s)); //  hat_s P 
  EC.ec_neg(tmp_ec);
  EC.ec_add(tilde_A_commit, tilde_A_commit, tmp_ec); // hat_sP - gamma Feld

  // Verif challenge

  tmp_mpz = H2(gg, pk, h, c, ECPointGroupCRefPair(A_commit, EC), tilde_h, tilde_V, ECPointGroupCRefPair(tilde_A_commit, EC)); 

  ret &= ( gamma == tmp_mpz) ;
 
  return ret;
}

void
CL_NIZK_RfSH_verify_multi(std::vector <bool> &res, std::vector <BICYCL::Mpz> gamma, std::vector <BICYCL::Mpz> hat_r, std::vector <BICYCL::Mpz> hat_s, std::vector <BICYCL::QFI> pk, std::vector <BICYCL::QFI> gg, std::vector <std::vector <BICYCL::QFI>> h, std::vector <std::vector <BICYCL::QFI>> c,  std::vector <BICYCL::ECPoint> &A_commit, std::vector <BICYCL::QFIPrecomp> gg_precomp, std::vector <BICYCL::QFIPrecomp> pk_precomp, BICYCL::CL_HSMqk CL, int start, int end) {


  for (int i = start; i < end ; i++) {

    res[i] = CL_NIZK_RfSH_verify(gamma[i], hat_r[i], hat_s[i], pk, gg, h[i], c[i], A_commit[i], gg_precomp, pk_precomp, CL);
    
  }
}

