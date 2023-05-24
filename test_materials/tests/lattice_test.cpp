/* lattice_test.cpp
*
* Unit tests for the Lattice used for LGMC
* All tests use catch2 unit test framework
*
*/

#include <gtest/gtest.h>
#include "../../LGMC/lattice.h"

// test that lattice has correct dimensions
TEST(LatticeIsConstructed, Lattice) {
  
   // create pointer to simple lattice
   Lattice *lattice = new Lattice(1, 3, 5, 7);
   ASSERT_EQ(lattice->xlo, 0);
   /*REQUIRE(lattice->xhi == 3);
   REQUIRE(lattice->ylo == 0);
   REQUIRE(lattice->yhi == 5);
   REQUIRE(lattice->zlo == 0);
   REQUIRE(lattice->zhi == 7);
   REQUIRE(lattice->latconst == 1);
   // x, y periodic and z non-periodic
   REQUIRE(lattice->sites.size() == 3*5*(7 + 1));
   REQUIRE(lattice->edges.size() == 3*5);*/


}


// test that each site has the correct neighbors



/*
// test add site
TEST_CASE("Site added to existing lattice", "[add_site]") {
   Lattice *lattice = new Lattice(1, 3, 5, 7);


   // z dimension grows
   lattice->add_site(2, 4, 8, 2, 4, 8, true, true, true);


   std::tuple<uint32_t, uint32_t, uint32_t> key = {2, 4, 8};
   int site_added = lattice->loc_map[key];


   REQUIRE(lattice->zhi == 8);
   REQUIRE(lattice->edges.size() == 3*5 + 1);
   REQUIRE(lattice->sites.size() == 3*5*(7+1) + 1);
  
   key = {2, 4, 7};
   int id = lattice->loc_map[key];
  
   REQUIRE(int(lattice->numneigh[site_added]) == 1);
   REQUIRE(int(lattice->idneigh[site_added][0]) == id);
}


// test delete site
TEST_CASE("Site deleted from existing lattice", "[delete_site]") {
   Lattice *lattice = new Lattice(1, 3, 5, 7);


   std::tuple<uint32_t, uint32_t, uint32_t> key = {2, 4, 7};
   int id = lattice->loc_map[key];


   lattice->delete_site(id);


   REQUIRE(lattice->zhi == 7);
   REQUIRE(lattice->edges.size() == 3*5 - 1);
   REQUIRE(lattice->sites.size() == 3*5*(7+1) - 1);
}


// test copy constructor
*/