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
}
