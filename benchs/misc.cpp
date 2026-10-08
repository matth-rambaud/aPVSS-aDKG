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

// precomp for faster exponentiation in basis g (adapted from BiCYCL for h in CL-HSM)
inline
void power_of_g (BICYCL::QFI &r, BICYCL::QFI g, BICYCL::Mpz x, BICYCL::QFIPrecomp g_precomp, BICYCL::CL_HSMqk CL)
  {
    CL.Cl_G().nupow (r, g, x, g_precomp);
  }


// check if an element is in the good group
bool
cg_check(BICYCL::QFI c, BICYCL::CL_HSMqk CL) {

  bool ret = true;
  
  ret &= c.discriminant() == CL.Cl_G().discriminant();
  ret &= CL.genus(c) == BICYCL::CL_HSMqk::Genus ({ 1, 1 });

  return ret;
}
