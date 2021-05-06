#include "gf2_polynomial_tests.h"
#include "gf2_polynomial.h"
#include "test_assert.h"

#include <sstream>

namespace {

void test_construction() {
  std::vector<uint8_t> c = {{1,0,0,1,1,0}};
  auto p = make_gf2_polynomial(c);
  TEST_EQ(p.coefficients.size(), 5);
}

void test_stream() {
  std::vector<uint8_t> c = {{1,0,0,1,1,0}};
  auto p = make_gf2_polynomial(c);
  std::stringstream ss;
  ss << p;
  TEST_EQ(ss.str(), std::string("X^4 + X^3 + 1"));
}

void test_stream_2() {
  gf2_polynomial g;
  g.coefficients.push_back(1);
  g.coefficients.push_back(2);
  g.coefficients.push_back(3);
  g.coefficients.push_back(4);
  g.coefficients.push_back(5);
  g.coefficients.push_back(6);
  std::stringstream ss;
  ss << g;
  TEST_EQ(ss.str(), std::string("X^4 + X^2 + 1"));
}

void test_degree() {
  gf2_polynomial g;
  g.coefficients.push_back(1);
  g.coefficients.push_back(2);
  g.coefficients.push_back(3);
  g.coefficients.push_back(4);
  g.coefficients.push_back(5);
  g.coefficients.push_back(6);
  TEST_EQ(degree(g), 4);
}

void test_equal() {
  gf2_polynomial g1;
  g1.coefficients.push_back(1);
  g1.coefficients.push_back(2);
  g1.coefficients.push_back(3);
  g1.coefficients.push_back(4);
  g1.coefficients.push_back(5);
  g1.coefficients.push_back(6);
  gf2_polynomial g2;
  g2.coefficients.push_back(1);
  g2.coefficients.push_back(0);
  g2.coefficients.push_back(1);
  g2.coefficients.push_back(0);
  g2.coefficients.push_back(1);
  TEST_ASSERT(g1==g2);
  
  gf2_polynomial g3;
  gf2_polynomial g4;
  g4.coefficients.push_back(2);
  g4.coefficients.push_back(4);
  g4.coefficients.push_back(6);
  g4.coefficients.push_back(8);
  TEST_ASSERT(g3==g4);
  TEST_ASSERT(g1!=g3);
  TEST_ASSERT(g1!=g4);
  TEST_ASSERT(g2!=g3);
  TEST_ASSERT(g2!=g4);
}

void test_add() {
  gf2_polynomial g1 = make_gf2_polynomial({{0,1,1,1}});
  gf2_polynomial g2 = make_gf2_polynomial({{1,0,1,0,1}});
  auto g3 = g1+g2;
  TEST_ASSERT(g3==make_gf2_polynomial({{1,1,0,1,1}}));
}

void test_sub() {
  gf2_polynomial g1 = make_gf2_polynomial({{0,1,1,1}});
  gf2_polynomial g2 = make_gf2_polynomial({{1,0,1,0,1}});
  auto g3 = g1-g2;
  TEST_ASSERT(g3==make_gf2_polynomial({{1,1,0,1,1}}));
}

void test_mul() {
  gf2_polynomial g1 = make_gf2_polynomial({{0,1,1,1}});
  gf2_polynomial g2 = make_gf2_polynomial({{1,0,1,0,1}});
  auto g3 = g1*g2;
  TEST_ASSERT(g3==make_gf2_polynomial({{0,1,1,0,1,0,1,1}}));
}

void test_derivative() {
  gf2_polynomial g1 = make_gf2_polynomial({{0,0,0,3,101,6,8}});
  gf2_polynomial g2 = derivative(g1);
  TEST_ASSERT(g2==make_gf2_polynomial({{0,0,1}}));
}

void test_xn() {
  gf2_polynomial g = make_xn(3);
  TEST_ASSERT(g==make_gf2_polynomial({{0,0,0,1}}));
  g = make_xn(2);
  TEST_ASSERT(g==make_gf2_polynomial({{0,0,1}}));
  g = make_xn(1);
  TEST_ASSERT(g==make_gf2_polynomial({{0,1}}));
  g = make_xn(0);
  TEST_ASSERT(g==make_gf2_polynomial({1}));
}

void test_euclidean_division() {
  gf2_polynomial a = make_gf2_polynomial({{0,0,0,1,1}});
  gf2_polynomial b = make_gf2_polynomial({{0,0,1}});
  auto div = euclidean_division(a, b);
  TEST_ASSERT(div.first==make_gf2_polynomial({{0,1,1}}));
  TEST_ASSERT(div.second==make_gf2_polynomial({0}));
  
  a = make_gf2_polynomial({{0,0,0,1,1,0,1,0}});
  b = make_gf2_polynomial({{0,1,1,1}});
  div = euclidean_division(a, b);
  //std::cout << a << " / " << b << " = (" << div.first << "," << div.second << ")\n";
  TEST_ASSERT(a == div.first*b + div.second);
  TEST_ASSERT(a/b == div.first);
  TEST_ASSERT(a%b == div.second);
  
  b = make_gf2_polynomial({{0,0,0,1,1,0,1,0}});
  a = make_gf2_polynomial({{0,1,1,1}});
  div = euclidean_division(a, b);
  //std::cout << a << " / " << b << " = (" << div.first << "," << div.second << ")\n";
  TEST_ASSERT(a == div.first*b + div.second);
  TEST_ASSERT(a/b == div.first);
  TEST_ASSERT(a%b == div.second);
}

void test_gcd() {
  gf2_polynomial a = make_gf2_polynomial({{0,0,0,1,1}});
  gf2_polynomial b = make_gf2_polynomial({{0,0,3,4}});
  gf2_polynomial g = gcd(a, b);
  TEST_ASSERT(g==make_gf2_polynomial({0,0,1}));
  TEST_ASSERT(a%g==make_gf2_polynomial({0}));
  TEST_ASSERT(b%g==make_gf2_polynomial({0}));
  
  a = make_gf2_polynomial({{0,0,0,1,1}});
  b = make_gf2_polynomial({{0,1,3,4}});
  g = gcd(a, b);
  TEST_ASSERT(g==make_gf2_polynomial({0,1,1}));
  TEST_ASSERT(a%g==make_gf2_polynomial({0}));
  TEST_ASSERT(b%g==make_gf2_polynomial({0}));
}

void test_hex_to_gf2_polynomial() {
  gf2_polynomial g = hex_to_gf2_polynomial("c");
  TEST_ASSERT(g==make_gf2_polynomial({0,0,1,1}));
  TEST_ASSERT(gf2_polynomial_to_hex(g)==std::string("c"));
  TEST_ASSERT(gf2_polynomial_to_hex(hex_to_gf2_polynomial("a466cfdc"))==std::string("a466cfdc"));
}

void test_sqrt() {
  gf2_polynomial g = make_gf2_polynomial({{1,0,1}});
  auto s = sqrt(g);
  TEST_ASSERT(s==make_gf2_polynomial({1,1}));
  g = hex_to_gf2_polynomial("a466cfdc");
  TEST_ASSERT(g == sqrt(power(g,2)));
}

void test_square_free_factorization() {
  gf2_polynomial g = hex_to_gf2_polynomial("c");
  auto R = square_free_factorization(g);
  auto p = R.front();
  for (int i = 1; i < R.size(); ++i)
  p = p * R[i];
  TEST_EQ(2, R.size());
  TEST_ASSERT(gf2_polynomial_to_hex(p) == gf2_polynomial_to_hex(g));
}

void test_distinct_degree_factorization() {
  gf2_polynomial g = make_gf2_polynomial({{0,0,1,1,0,1,0,0,1}});
  auto S = distinct_degree_factorization(g);
  auto p = make_xn(0);
  for (auto s : S) {
    //std::cout << s.first << std::endl;
    p = p*s.first;
  }
  TEST_ASSERT(p == g);
  TEST_EQ(3, S.size());
}

void test_equal_degree_factorization() {
  gf2_polynomial g = hex_to_gf2_polynomial("73af");
  uint64_t d = degree(g);
  auto factors = equal_degree_factorization(g, d/2);
  TEST_EQ(2, factors.size());
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[0]) == std::string("e5"));
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[1]) == std::string("83"));
  std::cout << gf2_polynomial_to_hex(factors[0]) << " " << gf2_polynomial_to_hex(factors[1]) << "\n";
}

void test_case_2() {
  gf2_polynomial g = hex_to_gf2_polynomial("738377c1");
  uint64_t d = degree(g);
  auto factors = equal_degree_factorization(g, d/2);
  //TEST_EQ(2, factors.size());
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[0]) == std::string("e5"));
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[1]) == std::string("83"));
  for (auto f : factors) {
    std::cout << gf2_polynomial_to_hex(f) << std::endl;  }
}

void test_case_3() {
  gf2_polynomial g = hex_to_gf2_polynomial("6677e20146508fb7");
  uint64_t d = degree(g);
  auto factors = equal_degree_factorization(g, d/2);
  //TEST_EQ(2, factors.size());
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[0]) == std::string("e5"));
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[1]) == std::string("83"));
  for (auto f : factors) {
    std::cout << gf2_polynomial_to_hex(f) << std::endl;  }
}

void test_case_4() {
//f3268b49 661859eb 0b324559 65ee6bda
  gf2_polynomial g = hex_to_gf2_polynomial("65ee6bda0b324559661859ebf3268b49");
  uint64_t d = degree(g);
  auto factors = equal_degree_factorization(g, d/2);
  //TEST_EQ(2, factors.size());
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[0]) == std::string("e5"));
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[1]) == std::string("83"));
  for (auto f : factors) {
    std::cout << gf2_polynomial_to_hex(f) << std::endl;  }
}

void test_case_5() {
//a91db473 fcea8db4 f3bb434a 8dba2f16 51abc87e 92c44759 5c1a16d3 6111c6f4
  gf2_polynomial g = hex_to_gf2_polynomial("6111c6f45c1a16d392c4475951abc87e8dba2f16f3bb434afcea8db4a91db473");
  uint64_t d = degree(g);
  auto factors = equal_degree_factorization(g, d/2);
  //TEST_EQ(2, factors.size());
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[0]) == std::string("e5"));
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[1]) == std::string("83"));
  for (auto f : factors) {
    std::cout << gf2_polynomial_to_hex(f) << std::endl;  }
}

void test_case_6() {
//4af6fc33 39029380 465c5267 c72f6a8b 0906e6d0 ca60550f 14a5e47c 42ad10fb 4a3bb446 bb74360a 5ea02b9c 23c68553 3fade253 e270ba24 39e141ad 6c38c43d
  gf2_polynomial g = hex_to_gf2_polynomial("6c38c43d39e141ade270ba243fade25323c685535ea02b9cbb74360a4a3bb44642ad10fb14a5e47cca60550f0906e6d0c72f6a8b465c5267390293804af6fc33");
  uint64_t d = degree(g);
  auto factors = equal_degree_factorization(g, d/2);
  //TEST_EQ(2, factors.size());
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[0]) == std::string("e5"));
  //TEST_ASSERT(gf2_polynomial_to_hex(factors[1]) == std::string("83"));
  for (auto f : factors) {
    std::cout << gf2_polynomial_to_hex(f) << std::endl;  }
}
} // namespace




void run_all_gf2_polynomial_tests() {
  test_construction();
  test_stream();
  test_stream_2();
  test_degree();
  test_equal();
  test_add();
  test_sub();
  test_mul();
  test_derivative();
  test_xn();
  test_euclidean_division();
  test_gcd();
  test_hex_to_gf2_polynomial();
  test_sqrt();
  test_square_free_factorization();
  test_distinct_degree_factorization();
  test_equal_degree_factorization();
  test_case_2();
  test_case_3();
  test_case_4();
  test_case_5();
  test_case_6();
}
