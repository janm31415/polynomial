#ifndef GF2_POLYNOMIAL_H
#define GF2_POLYNOMIAL_H

#include <iostream>
#include <vector>
#include <cmath>

struct gf2_polynomial {
  std::vector<uint8_t> coefficients;
};

std::vector<uint8_t> simplify_gf2_coefficients(const std::vector<uint8_t>& coeff) {
  std::vector<uint8_t> c(coeff);
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

inline gf2_polynomial make_gf2_polynomial(const std::vector<uint8_t>& coef) {
  gf2_polynomial g;
  g.coefficients = simplify_gf2_coefficients(coef);
  return g;
}

inline gf2_polynomial hex_to_gf2_polynomial(std::string hexadecimal_number) {
  std::vector<uint8_t> coeff;
  coeff.reserve(4*hexadecimal_number.length());
  while (!hexadecimal_number.empty()) {
    char ch = hexadecimal_number.back();
    int i = 0;
    if (ch >= '0' && ch <= '9')
      i = (int)(ch-'0');
    else if (ch >= 'A' && ch <= 'F')
      i = (int)(ch-'A'+10);
    else if (ch >= 'a' && ch <= 'f')
      i = (int)(ch-'a'+10);
    else throw std::runtime_error("make_gf2_polynomial: input string is not a hexadecimal number!");
    coeff.push_back(i&1?1:0);
    coeff.push_back(i&2?1:0);
    coeff.push_back(i&4?1:0);
    coeff.push_back(i&8?1:0);
    hexadecimal_number.pop_back();
  }
  return make_gf2_polynomial(coeff);
}

inline std::string gf2_polynomial_to_hex(const gf2_polynomial& g) {
  std::string s;
  const size_t sz = g.coefficients.size();
  s.reserve(sz/4+1);
  for (size_t i = 0; i < sz; i+=4) {
    int i0 = g.coefficients[i];
    int i1 = i+1 < sz ? g.coefficients[i+1]:0;
    int i2 = i+2 < sz ? g.coefficients[i+2]:0;
    int i3 = i+3 < sz ? g.coefficients[i+3]:0;
    int h = i0+2*i1+4*i2+8*i3;
    if (h<10)
      s.push_back(h+'0');
    else
      s.push_back(h-10+'a');
  }
  std::reverse(s.begin(), s.end());
  return s;
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
  std::vector<uint8_t> coeff;
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
  std::vector<uint8_t> coeff;
  coeff.reserve(p.coefficients.size());
  for (size_t i = 1; i < p.coefficients.size(); ++i) {
    if ((i&1)==0 || (p.coefficients[i]&1)==0)
      coeff.push_back(0);
    else
      coeff.push_back(1);
  }
  return make_gf2_polynomial(coeff);
}

inline void minus_b_times_xn(gf2_polynomial& r, const gf2_polynomial& b, uint64_t n)
  {
  if (r.coefficients.size() < b.coefficients.size()+n)
    r.coefficients.resize(b.coefficients.size()+n, 0);
  for (uint64_t i = 0; i < b.coefficients.size(); ++i) {
    r.coefficients[i+n] += b.coefficients[i];
    r.coefficients[i+n] &= 1;
  }
  while (!r.coefficients.empty() && ((r.coefficients.back()&1)==0))
    r.coefficients.pop_back();
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
  q.coefficients.resize(qsize,0);
  uint64_t deg_r = degree(r);
  while(!r.coefficients.empty() && deg_r >= d) {
    uint64_t deg_b = degree(b);
    uint64_t n = deg_r - deg_b;
    q.coefficients[n] = 1;
    minus_b_times_xn(r, b, n);
    deg_r = degree(r);
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
    std::swap(a.coefficients, b.coefficients);
    std::swap(r.coefficients, b.coefficients);
    //a = b;
    //b = r;
    r = a % b;
  }
  return b;
}

inline gf2_polynomial power(const gf2_polynomial& a, int p) {
  gf2_polynomial result = make_gf2_polynomial({1});
  for (int i = 1; i <= p; ++i)
    result = result * a;
  return result;
}

inline gf2_polynomial sqrt(const gf2_polynomial& a) {
  gf2_polynomial result;
  for (int i = 0; i < a.coefficients.size(); ++i)
    if (i%2==0) {
      result.coefficients.push_back(a.coefficients[i]);
    }
  return result;
}

inline gf2_polynomial make_random_gf2_polynomial(uint64_t n) {
  gf2_polynomial p;
  for (uint64_t i = 0; i <= n; ++i)
    p.coefficients.push_back(rand()&1);
  return simplify(p);
}

/*
inline std::pair<gf2_polynomial, gf2_polynomial> split_in_square_free_part(const gf2_polynomial& a) {
  gf2_polynomial d = derivative(a);
  gf2_polynomial g = gcd(a, d);
  gf2_polynomial r = a/g;
  return std::pair<gf2_polynomial, gf2_polynomial>(r, g);
}
*/

//source: https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
inline std::vector<gf2_polynomial> square_free_factorization(const gf2_polynomial& f) {
  std::vector<gf2_polynomial> R;
  R.push_back(make_gf2_polynomial({1}));
  
  //Make w be the product (without multiplicity) of all factors of f that have
  //multiplicity not divisible by p
    
  auto c = gcd(f, derivative(f));
  auto w = f/c;
  
  // Step 1: Identify all factors in w
  int i = 1;
  while (w != R.front()) {
    auto y = gcd(w, c);
    auto fac = w/y;
    R.push_back(power(fac, i));
    w = y;
    c = c/y;
    ++i;
  }
  // c is now the product (with multiplicity) of the remaining factors of f
   
  // Step 2: Identify all remaining factors using recursion
  // Note that these are the factors of f that have multiplicity divisible by p
  if (c != R.front()) {
    c = sqrt(c);
    auto R2 = square_free_factorization(c);
    for (const auto& factor : R2) {
      if (factor != R.front())
        R.push_back(power(factor, 2));
      }
  }
  if (R.size()>1)
    R.erase(R.begin());
  return R;
}

//source: https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
std::vector<std::pair<gf2_polynomial, uint64_t>> distinct_degree_factorization(const gf2_polynomial& f) {
/*
    Input: A monic square-free polynomial f in GF2
    Output: The set of all pairs (g, d), such that
             f has an irreducible factor of degree d and
             g is the product of all monic irreducible factors of f of degree d.
*/
  std::vector<std::pair<gf2_polynomial, uint64_t>> S;
  int i = 1;
  auto fstar = f;
  auto unit = make_xn(0);
  while (degree(fstar)>=2*i) {
    auto g = gcd(fstar, make_xn((uint64_t)pow(2,i)) - make_xn(1));
    if (g != unit) {
      S.emplace_back(g, i);
      fstar = fstar/g;
    }
    ++i;
  }
  if (fstar != unit) {
    S.emplace_back(fstar, degree(fstar));
  }
  if (S.empty())
    S.emplace_back(f, 1);
  return S;
}

std::vector<gf2_polynomial> equal_degree_factorization(const gf2_polynomial& f, uint64_t d) {
/*
Input: A monic square free polynomial f in GF2 of degree n = rd, which
       has r >= 2 irreducible factors each of degree d.
Output: The set of monic irreducible factors of f.

Source for p=2: https://math.stackexchange.com/questions/1636518/how-do-i-apply-the-cantor-zassenhaus-algorithm-to-mathbbf-2
*/
  std::vector<gf2_polynomial> factors;
  factors.push_back(f);
  auto unit = make_xn(0);
  
  uint64_t n = degree(f);
  
  uint64_t r = n/d;
  
  while (factors.size() < r) {
    auto h = make_random_gf2_polynomial(n-1);
    auto g = h;
    //g = h + h^2 + h^4 + ... + h^(2^(d-1))
    auto last_term = h;
    for (int j = 1; j < d; ++j) {
      last_term = (last_term*last_term) % f;
      if (last_term.coefficients.empty())
        break;
      g = g + last_term;
      }
    //g = g%f;
    if (g.coefficients.empty())
      continue;
    for (size_t i = 0; i < factors.size(); ++i) {
      const auto& u = factors[i];
      if (degree(u)>d) {
        auto gcd_g_u = gcd(g, u);
        if (gcd_g_u != unit && gcd_g_u != u) {
          factors.push_back(u/gcd_g_u);
          factors[i] = gcd_g_u;  
        }
      }
    }
  }
  
  return factors;
}

#endif
