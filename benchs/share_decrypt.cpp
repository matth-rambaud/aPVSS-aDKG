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

// encryption subfunction for threads: compute element start to end of ciphertext
void
partial_encryption(std::vector <BICYCL::QFI> &h_ciphertext, std::vector <BICYCL::QFI> &c_ciphertext, std::vector <BICYCL::QFI> gg, std::vector <BICYCL::QFI> pk, std::vector <BICYCL::Mpz> s, BICYCL::Mpz r_ciphertext, std::vector <BICYCL::QFIPrecomp> gg_precomp, std::vector <BICYCL::QFIPrecomp> pk_precomp, BICYCL::CL_HSMqk CL, unsigned int start, unsigned int end)
{
  BICYCL::QFI tmp_qfi;
  
  for (unsigned int i = start; i < end; i++) {
    if (i <= K) {
      // g0^r,...,gK^r
      power_of_g(h_ciphertext[i], gg[i], r_ciphertext, gg_precomp[i], CL);  
    }

    else {
      // f^s1 pk0^r,...f^sn pk_{n-1}^r
      power_of_g(c_ciphertext[i-K-1], pk[i-K-1], r_ciphertext, pk_precomp[i-K-1], CL);
    tmp_qfi = CL.power_of_f(s[i-K]);
    CL.Cl_G().nucomp(c_ciphertext[i-K-1], c_ciphertext[i-K-1], tmp_qfi);
    }
  }
}



// share function
// input s \in Fq
// output (a0, ..., at) coeff Eval Poly
// share (s1, ... ,sn) shares of s
// (h,c) ciphertext
// r randomness ciphertext

void
share(std::vector <BICYCL::Mpz> &shares, std::vector <BICYCL::Mpz> &coeffs, BICYCL::Mpz &r_ciphertext,  std::vector <BICYCL::QFI> &h_ciphertext, std::vector <BICYCL::QFI> &c_ciphertext, BICYCL::Mpz secret, std::vector <BICYCL::QFI> gg, std::vector <BICYCL::QFI> pk, BICYCL::Mpz bound_expo, std::vector <BICYCL::QFIPrecomp> gg_precomp, std::vector <BICYCL::QFIPrecomp> pk_precomp, BICYCL::CL_HSMqk CL,BICYCL::RandGen &randgen)
{

  coeffs[0] = secret;
  
  for (unsigned int i = 1; i < t_param + 1 ; i++) {
    coeffs[i] = randgen.random_mpz (CL.q());
  }
  
  // Compute evalutations with Horner's method in i=1,...,n
  // shares[i] = a0 + a1 (i) + ··· + at (i)^t
  //=    a0 +  (        ...  + (at-1   + at i)*i ... )*i
  
    
  for (unsigned int i = 1; i < n + 1; i++) {
    shares[i] = coeffs[t_param];
    for (int k = t_param - 1; k >= 0; k--)
      {
	BICYCL::Mpz::mul(shares[i], shares[i], i);
	BICYCL::Mpz::add(shares[i], shares[i], coeffs[k]);
	BICYCL::Mpz::mod(shares[i], shares[i], CL.q());
      }
  }

  // Encryption

  r_ciphertext = randgen.random_mpz (bound_expo);
  
  /* without threads 
  // g0^r,...,gK^r
  
  for (unsigned int j = 0; j < K + 1; j++) {
    power_of_g(h_ciphertext[j], gg[j], r_ciphertext, gg_precomp[j], CL);  
  }

  // f^s1 pk0^r,...f^sn pk_{n-1}^r
  
  for (unsigned int i = 0; i < n; i++) {
    power_of_g(c_ciphertext[i], pk[i], r_ciphertext, pk_precomp[i], CL);
    tmp_qfi = CL.power_of_f(s[i+1]);
    CL.Cl_G().nucomp(c_ciphertext[i], c_ciphertext[i], tmp_qfi);
  }
  */

  // with threads
  auto th = new std::thread[nb_threads];
  unsigned int nb = K + n + 1; 
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

  first = 0;
  for (unsigned int tt = 0; tt < local_nb_threads - 1; tt++) {
    last = first + quotient + (tt < remainder ? 1 : 0);
    th[tt] = std::thread( partial_encryption, std::ref(h_ciphertext), std::ref(c_ciphertext), gg, pk, shares, r_ciphertext, gg_precomp, pk_precomp, CL, first, last);
    first = last;
  }

  partial_encryption(h_ciphertext, c_ciphertext, gg, pk, shares, r_ciphertext, gg_precomp, pk_precomp, CL, first, nb);
  
  for (unsigned int tt = 0; tt < local_nb_threads-1; tt++) {
      th[tt].join();
  }
}

void
decrypt(BICYCL::QFI &A, BICYCL::Mpz &plaintext, std::vector <BICYCL::QFI> h_ciphertext, BICYCL::QFI c_ciphertext_i, std::vector <BICYCL::Mpz> sk, BICYCL::CL_HSMqk CL)
{
  BICYCL::QFI tmp_qfi;

  CL.Cl_G().nupow(A, h_ciphertext, sk, 6); // multi-expo \prod h_i^ski
  CL.Cl_G().nucompinv(tmp_qfi, c_ciphertext_i, A);
  plaintext = CL.dlog_in_F(tmp_qfi);
}

