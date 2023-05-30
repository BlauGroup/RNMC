/* lattice_test.cpp
*
* Unit tests for the Lattice used for LGMC
* All tests use catch2 unit test framework
*
*/

#include <gtest/gtest.h>
#include "../../LGMC/lattice.h"

// test that lattice has correct dimensions
TEST(create_basic_lattice, Lattice) {
  
   // create pointer to simple lattice
   Lattice *lattice = new Lattice(1, 3, 5, 7);
   ASSERT_EQ(lattice->xlo, 0);
   ASSERT_EQ(lattice->xhi, 3);
   ASSERT_EQ(lattice->ylo, 0);
   ASSERT_EQ(lattice->yhi, 5);
   ASSERT_EQ(lattice->zlo, 0);
   ASSERT_EQ(lattice->zhi, 7);
   ASSERT_EQ(lattice->latconst, 1);
   
   // x, y periodic and z non-periodic
   ASSERT_EQ(static_cast<int>(lattice->sites.size()), 3*5*(7 + 1));
   ASSERT_EQ(static_cast<int>(lattice->edges.size()), 3*5);


}


// test that each site has the correct neighbors


// test add site
TEST(add_site, add_site) {
   Lattice *lattice = new Lattice(1, 3, 5, 7);

   // z dimension grows
   lattice->add_site(2, 4, 8, 2, 4, 8, true, true, true);

   std::tuple<uint32_t, uint32_t, uint32_t> key = {2, 4, 8};
   int site_added = lattice->loc_map[key];


   ASSERT_EQ(lattice->zhi, 8);
   ASSERT_EQ(static_cast<int>(lattice->edges.size()), 3*5 + 1);
   ASSERT_EQ(static_cast<int>(lattice->sites.size()), 3*5*(7+1) + 1);
  
   key = {2, 4, 7};
   int id = lattice->loc_map[key];
  
   ASSERT_EQ(int(lattice->numneigh[site_added]), 1);
   ASSERT_EQ(int(lattice->idneigh[site_added][0]), id);
}


// test delete site
TEST(delete_site, delete_site) {
   Lattice *lattice = new Lattice(1, 3, 5, 7);

   std::tuple<uint32_t, uint32_t, uint32_t> key = {2, 4, 7};
   int id = lattice->loc_map[key];

   lattice->delete_site(id);

   ASSERT_EQ(lattice->zhi, 7);
   ASSERT_EQ(static_cast<int>(lattice->edges.size()), 3*5 - 1);
   ASSERT_EQ(static_cast<int>(lattice->sites.size()), 3*5*(7+1) - 1);
}


// test copy constructor
