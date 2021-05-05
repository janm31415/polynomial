#ifndef GF2_POLYNOMIAL_H
#define GF2_POLYNOMIAL_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <bitset>
#include <cmath>
#include <stdint.h>

struct gf2_polynomial {
  std::vector<uint64_t> coefficients;
};

std::vector<uint64_t> simplify_gf2_coefficients(const std::vector<uint64_t>& coeff) {
  std::vector<uint64_t> c(coeff);
  while (!c.empty() && ((c.back()&1)==0))
    c.pop_back();
  for (auto& v : c)
    v = v&1;
  return c;
}

inline gf2_polynomial make_gf2_polynomial(const std::vector<uint64_t>& coef) {
  gf2_polynomial g;
  g.coefficients = simplify_gf2_coefficients(coef);
  return g;
}

std::ostream& operator<<(std::ostream& s, const gf2_polynomial& p) {
  bool first = true;
  if (p.coefficients.empty())
    s << "0";
  int count = p.coefficients.size()-1;
  for (auto rit = p.coefficients.rbegin(); rit != p.coefficients.rend(); ++rit,--count) {
    if(*rit & 1) {
      if (!first)
        s << " + ";
      if (count)
        s << "X^"<<count;
      else
        s << "1";
      first = false;
    }
  }
  return s;
}

inline uint64_t degree(const gf2_polynomial& p) {
  if (p.coefficients.empty())
    return 0;
  uint64_t deg = p.coefficients.size()-1;
  while (deg && ((p.coefficients[deg]&1)==0))
    --deg;
  return deg;
}

#endif
