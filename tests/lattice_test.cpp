/* 
*
* Unit tests for the Lattice used for LGMC
* All tests use googletest unit test framework
*
*/

#include <gtest/gtest.h>
#include "../LGMC/lattice.h"

TEST(lattice_test, InitalizationWorks) {
  
   // create pointer to simple lattice
   Lattice *lattice = Lattice(1, 3, 5, 7);
   EXPECT_EQ(lattice.xlo, 0);
   EXPECT_EQ(lattice.xhi, 3);
   EXPECT_EQ(lattice.ylo, 0);
   EXPECT_EQ(lattice.yhi, 5);
   EXPECT_EQ(lattice.zlo, 0);
   EXPECT_EQ(lattice.zhi, 7);
   EXPECT_EQ(lattice.latconst, 1);
   
   // x, y periodic and z non-periodic
   EXPECT_EQ(static_cast<int>(lattice.sites.size()), 3*5*(7 + 1));
   EXPECT_EQ(static_cast<int>(lattice.edges.size()), 3*5);

   delete lattice;
}


// test that each site has the correct neighbors

// test add site
TEST(lattice_test, AddSite) {
   Lattice *lattice = Lattice(1, 3, 5, 7);

   // z dimension grows
   lattice.add_site(2, 4, 8, true, true, true);

   std::tuple<uint32_t, uint32_t, uint32_t> key = {2, 4, 8};
   int site_added = lattice.loc_map[key];


   ASSERT_EQ(lattice.zhi, 8);
   ASSERT_EQ(static_cast<int>(lattice.edges.size()), 3*5 + 1);
   ASSERT_EQ(static_cast<int>(lattice.sites.size()), 3*5*(7+1) + 1);
  
   key = {2, 4, 7};
   int id = lattice.loc_map[key];
  
   ASSERT_EQ(int(lattice.numneigh[site_added]), 1);
   ASSERT_EQ(int(lattice.idneigh[site_added][0]), id);

   delete lattice;
}


// test delete site
TEST(lattice_test, DeleteSite) {
   Lattice *lattice = Lattice(1, 3, 5, 7);

   std::tuple<uint32_t, uint32_t, uint32_t> key = {2, 4, 7};
   int id = lattice.loc_map[key];

   lattice.delete_site(id);

   ASSERT_EQ(lattice.zhi, 7);
   ASSERT_EQ(static_cast<int>(lattice.edges.size()), 3*5 - 1);
   ASSERT_EQ(static_cast<int>(lattice.sites.size()), 3*5*(7+1) - 1);

   delete lattice;
}

// test copy constructor
TEST(lattice_test, CopyConstructor) {
   Lattice *lattice = Lattice(1, 3, 5, 7);

   Lattice *new_lattice = nullptr;

   new_lattice = lattice;
   lattice = nullptr;

   EXPECT_EQ(new_lattice.xlo, 0);
   EXPECT_EQ(new_lattice.xhi, 3);
   EXPECT_EQ(new_lattice.ylo, 0);
   EXPECT_EQ(new_lattice.yhi, 5);
   EXPECT_EQ(new_lattice.zlo, 0);
   EXPECT_EQ(new_lattice.zhi, 7);
   EXPECT_EQ(new_lattice.latconst, 1);
   
   // x, y periodic and z non-periodic
   EXPECT_EQ(static_cast<int>(new_lattice.sites.size()), 3*5*(7 + 1));
   EXPECT_EQ(static_cast<int>(new_lattice.edges.size()), 3*5);


   delete new_lattice;
}
