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

// w-NAF of k  Algo 3.35  guide to ECC. pad to get size pad
inline void WNAF(std::vector <long> &W,  BICYCL::Mpz k, int w, size_t pad) {

  BICYCL::Mpz aux;

  W = {};
  
  while (k >= 1L) { 
    if (k.tstbit (0)) {   // k%2
      if (w == 1) {
	BICYCL::Mpz::mod2k(aux, k, w); // because centered gives -1, bug bicycl?
      }
      else {
	BICYCL::Mpz::mod2k_centered(aux, k, w);  // aux = k mod_c 2^w
      }
	      
      W.push_back((long) aux);
      BICYCL::Mpz::sub(k, k, aux);   // k = k - aux
    }
    else {
      W.push_back(0L);
    }
    BICYCL::Mpz::divby2(k, k);  // k = k//2
  }

  // pad to get size = pad
  while (W.size() < pad) {
      W.push_back(0L);
  }
}



// Multi-expo EC (deprecated in OpenSSL)
// C = sum_0^{nb-1} expos[i] gens[i] with WNAF. Assume expos of bitsize q_size and C = 0_EC
void
EC_multi_exp(BICYCL::ECPoint &C, std::vector <BICYCL::ECPoint> &gens,  std::vector <BICYCL::Mpz> &expos, unsigned int nb,  const BICYCL::ECGroup &EC) {
  auto W = new std::vector <long>[nb];

  // compute the NAFs
  for (unsigned int i = 0; i < nb; i++) {
    WNAF(W[i], expos[i], w_ecmulti_exp, q_size + 1); 
  }

  std::vector <std::vector <BICYCL::ECPoint>> P ;
  //size [nb][nb_precomp_w_ecmulti_exp];
  BICYCL::ECPoint aux (EC);
  BICYCL::ECPoint aux2 (EC);

  for (unsigned int i = 0; i < nb; i++) {
    P.push_back (std::vector<BICYCL::ECPoint>());
    P[i].push_back(BICYCL::ECPoint(EC, gens[i])); // P[i][0] = gens[i]
    EC.ec_add(aux, P[i][0], P[i][0]);   // aux =2 P
    for (unsigned int j = 1; j < nb_precomp_w_ecmulti_exp; j++) {
      EC.ec_add(aux2, P[i][j-1], aux); // aux 2 =(2(j-1) + 1) P_i + 2 P_i = (2j+1) P_i
      P[i].push_back(BICYCL::ECPoint(EC, aux2)); //P[i][j] = aux2
    }
  }
  
  for (int i = q_size; i > -1; i--) {
    EC.ec_add(C, C, C);  // C <- 2C
    for (unsigned int j = 0; j < nb; j++) {
      if (W[j][i]) {
	aux = P[j][(abs(W[j][i])-1)/2];
	if (W[j][i] < 0){
	  EC.ec_neg(aux);
	}
	EC.ec_add(C, C, aux);
      } // C <- C +-  P[j][(abs(W[j][i])-1)/2]
    }
  }
}
