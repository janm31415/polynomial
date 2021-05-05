#include "gf2_polynomial_tests.h"
#include "gf2_polynomial.h"
#include "test_assert.h"

#include <sstream>

namespace {

  void test_construction() {
    std::vector<uint64_t> c = {{1,0,0,1,1,0}};
    auto p = make_gf2_polynomial(c);
    TEST_EQ(p.coefficients.size(), 5);
  }

  void test_stream() {
    std::vector<uint64_t> c = {{1,0,0,1,1,0}};
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
  
  //46508fb7  6677e201
  void test_factorization() {
  /*
    auto product = hex_to_gf2_polynomial("b0c152f9")*hex_to_gf2_polynomial("ebf2831f");
    gf2_polynomial g = hex_to_gf2_polynomial("6677e20146508fb7");
    std::cout << "Product: " << gf2_polynomial_to_hex(product) << std::endl;
    std::cout << "Input: " << gf2_polynomial_to_hex(g) << std::endl;
    auto R = square_free_factorization(g);
    for (auto f : R) {
      std::cout << gf2_polynomial_to_hex(f) << std::endl;
    }
    auto p = R.front();
    for (int i = 1; i < R.size(); ++i)
      p = p * R[i];
    std::cout << "Product of factors: " << gf2_polynomial_to_hex(p) << std::endl;
    std::cout << "Original polynomial: " << gf2_polynomial_to_hex(g) << std::endl;
    */
    //gf2_polynomial g = hex_to_gf2_polynomial("6677e20146508fb7");
    //gf2_polynomial s = sqrt(g);
    //std::cout << "Original polynomial: " << gf2_polynomial_to_hex(g) << std::endl;
    //std::cout << "sqrt^2: " << gf2_polynomial_to_hex(power(s,2)) << std::endl;
    gf2_polynomial g = hex_to_gf2_polynomial("6677e20146508fb7");
    //gf2_polynomial g = make_gf2_polynomial({{0,0,1,1,0,1,0,0,1}});
    //auto S = distinct_degree_factorization(g);
    //auto p = make_xn(0);
    //for (auto s : S) {
    //  std::cout << gf2_polynomial_to_hex(s.first) << std::endl;
    //  p = p*s.first;
    //}
    //std::cout << p << std::endl;
  }

}




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
  test_factorization();
}
