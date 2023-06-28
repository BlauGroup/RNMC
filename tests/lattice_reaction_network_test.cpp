/* -------------------------------------------------------------------------
*
* Unit tests for lattice reaction network
* All tests use catch2 unit test framework
*
------------------------------------------------------------------------- */

#include "../core/sql.h"
#include "../LGMC/lattice_reaction_network.h"

#include <gtest/gtest.h>

#include <string>

TEST(lattice_reaction_network_test, InitialState) {

   // static lattice
   std::string model_database_file = "../examples/LGMC/CO_oxidation/rn.sqlite";
   std::string initial_state_database_file = "../examples/LGMC/CO_oxidation/initial_state.sqlite";

   SqlConnection model_database = SqlConnection(model_database_file,
                              SQLITE_OPEN_READWRITE);
   SqlConnection initial_state_database = SqlConnection(initial_state_database_file,
                              SQLITE_OPEN_READWRITE);

   LatticeParameters parameters = {.latconst = 1,
                           .boxxhi = 50,
                           .boxyhi = 50,
                           .boxzhi = 2,
                           .temperature = 300,
                           .g_e = -0.5, .is_add_sites = false,
                           .charge_transfer_style = ChargeTransferStyle::BUTLER_VOLMER}; 


   LatticeReactionNetwork static_LGMC_ = LatticeReactionNetwork(model_database,
                                       initial_state_database,
                                       parameters);

   // dynamic lattice 
   model_database_file = "../examples/LGMC/SEI/rn.sqlite";
   initial_state_database_file = "../examples/LGMC/SEI/initial_state.sqlite";

   model_database = SqlConnection(model_database_file,
                              SQLITE_OPEN_READWRITE);
   initial_state_database = SqlConnection(initial_state_database_file,
                              SQLITE_OPEN_READWRITE);

   parameters = {.latconst = 1,
               .boxxhi = 100,
               .boxyhi = 100,
               .boxzhi = 2,
               .temperature = 300,
               .g_e = -2.1, .is_add_sites = true,
               .charge_transfer_style = ChargeTransferStyle::MARCUS}; 


   LatticeReactionNetwork dynamic_LGMC_ = LatticeReactionNetwork(model_database,
                                          initial_state_database,
                                          parameters);
   
   // check static lattice
   EXPECT_EQ(static_LGMC_.initial_state[0], 0);
   EXPECT_EQ(static_LGMC_.initial_state[1], 2500);
   EXPECT_EQ(static_LGMC_.initial_state[2], 0);
   EXPECT_EQ(static_LGMC_.initial_state[3], 0);
   EXPECT_EQ(static_LGMC_.initial_state[4], 15000);

   // check dynamic lattice
   for(int i = 0; i < static_cast<int>(dynamic_LGMC_.initial_state.size()); i++) {
      if(i == 3) {
         EXPECT_EQ(dynamic_LGMC_.initial_state[3], 10000);
      }
      else {
         EXPECT_EQ(dynamic_LGMC_.initial_state[i], 0);
      }
   }

}






