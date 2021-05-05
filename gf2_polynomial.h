#ifndef GF2_POLYNOMIAL_H
#define GF2_POLYNOMIAL_H

#include <iostream>
#include <vector>

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

inline gf2_polynomial simplify(const gf2_polynomial& g) {
  gf2_polynomial ret;
  ret.coefficients = simplify_gf2_coefficients(g.coefficients);
  return ret;
}

inline gf2_polynomial make_gf2_polynomial(const std::vector<uint64_t>& coef) {
  gf2_polynomial g;
  g.coefficients = simplify_gf2_coefficients(coef);
  return g;
}

inline gf2_polynomial make_xn(uint64_t n) {
  gf2_polynomial g;
  for (size_t i = 0; i < n; ++i)
    g.coefficients.push_back(0);
  g.coefficients.push_back(1);
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

inline bool operator == (const gf2_polynomial& a, const gf2_polynomial& b) {
  uint64_t d1 = degree(a);
  uint64_t d2 = degree(b);
  if (d1 != d2)
    return false;
  if (a.coefficients.empty())
    return b.coefficients.empty() || (b.coefficients.front()&1)==0;
  if (b.coefficients.empty())
    return (a.coefficients.front()&1)==0;
  for (uint64_t i = 0; i <= d1; ++i) {
    if ((a.coefficients[i]&1) != (b.coefficients[i]&1))
      return false;
  }
  return true;
}

inline bool operator != (const gf2_polynomial& a, const gf2_polynomial& b) {
  return !(a==b);
}

inline gf2_polynomial operator + (const gf2_polynomial& a, const gf2_polynomial& b) {
  if (a.coefficients.size()>=b.coefficients.size()) {
    auto coef = a.coefficients;
    for (size_t i=0; i < b.coefficients.size(); ++i)
      coef[i] += b.coefficients[i];
    return make_gf2_polynomial(coef);
  } else {
    auto coef = b.coefficients;
    for (size_t i=0; i < a.coefficients.size(); ++i)
      coef[i] += a.coefficients[i];
    return make_gf2_polynomial(coef);
  }
}

inline gf2_polynomial operator - (const gf2_polynomial& a, const gf2_polynomial& b) {
  return a+b;
}


/*
    carry-less multiplication from http://bitmath.blogspot.com/2013/05/carryless-multiplicative-inverse.html
    uint r = 0;
    while (b != 0)
    {
        if ((a & 1) != 0)
            r ^= b;      // carryless addition is xor
        a >>= 1;
        b <<= 1;
    }
    return r;
*/
inline gf2_polynomial operator * (const gf2_polynomial& a, const gf2_polynomial& b) {
  uint64_t a_index = 0;
  std::vector<uint64_t> coeff;
  coeff.reserve(a.coefficients.size()*b.coefficients.size());
  while (a_index < a.coefficients.size()) {
    if (a.coefficients[a_index]&1) {
      coeff.resize(b.coefficients.size()+a_index);
      for (size_t i = 0; i < b.coefficients.size(); ++i)
        coeff[i+a_index] ^= b.coefficients[i];
    }
    ++a_index;
  }
  return make_gf2_polynomial(coeff);
}

inline gf2_polynomial derivative(const gf2_polynomial& p) {
  std::vector<uint64_t> coeff;
  coeff.reserve(p.coefficients.size());
  for (size_t i = 1; i < p.coefficients.size(); ++i) {
    if ((i&1)==0 || (p.coefficients[i]&1)==0)
      coeff.push_back(0);
    else
      coeff.push_back(1);
  }
  return make_gf2_polynomial(coeff);
}

/*
The Euclidean division provides two polynomials q(x), the quotient and r(x), the remainder such that
a(x)=q0(x)b(x)+r0(x) and deg⁡(r0(x)) < deg⁡(b(x))
*/
inline std::pair<gf2_polynomial, gf2_polynomial> euclidean_division(const gf2_polynomial& a, const gf2_polynomial& b) {
  gf2_polynomial r(a);
  gf2_polynomial q;
  auto d = degree(b);
  auto deg_a = degree(a);
  uint64_t qsize = deg_a>d?deg_a-d+1:1;
  q.coefficients.resize(qsize);
  while(degree(r) >= d) {
    uint64_t n = degree(r) - degree(b);
    q.coefficients[n] = (r.coefficients[degree(r)] / b.coefficients[degree(b)])&1;
    if (q.coefficients[n])
      r = r - b*make_xn(n);
  }
  return std::make_pair(simplify(q), simplify(r));
}

inline gf2_polynomial operator / (const gf2_polynomial& a, const gf2_polynomial& b) {
  return euclidean_division(a,b).first;
}

inline gf2_polynomial operator % (const gf2_polynomial& a, const gf2_polynomial& b) {
  return euclidean_division(a,b).second;
}

inline gf2_polynomial gcd(gf2_polynomial a, gf2_polynomial b) {
  if (degree(a)<degree(b))
    std::swap(a.coefficients, b.coefficients);
  gf2_polynomial r = a % b;
  while(!r.coefficients.empty()) {
    a = b;
    b = r;
    r = a % b;
  }
  return b;
}

#endif
