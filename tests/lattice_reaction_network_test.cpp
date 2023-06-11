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

class LatticeReactionNetworkTest : public ::testing::Test {
   protected:
   void SetUp() override {
      
      // static lattice
      std::string model_database_file = "../examples/LGMC/CO_oxidation(static)/network.sqlite";
      std::string initial_state_database_file = "../examples/LGMC/CO_oxidation(static)/state.sqlite";

      SqlConnection model_database = SqlConnection(model_database_file,
                                 SQLITE_OPEN_READWRITE);
      SqlConnection initial_state_database = SqlConnection(initial_state_database_file,
                                 SQLITE_OPEN_READWRITE);

      LatticeParameters parameters{.latconst = 1,
                              .boxxhi = 50,
                              .boxyhi = 50,
                              .boxzhi = 2,
                              .temperature = 300,
                              .g_e = -0.5, .is_add_sites = false,
                              .charge_transfer_style = ChargeTransferStyle::BUTLER_VOLMER,
                              .lattice_fill = "none"}; 


      static_LGMC_ = LatticeReactionNetwork(model_database,
                                          initial_state_database,
                                          parameters);

      // dynamic lattice 
      std::string model_database_file = "../examples/LGMC/SEI(dynamic)/network.sqlite";
      std::string initial_state_database_file = "../examples/LGMC/SEI(dynamic)/state.sqlite";

      SqlConnection model_database = SqlConnection(model_database_file,
                                 SQLITE_OPEN_READWRITE);
      SqlConnection initial_state_database = SqlConnection(initial_state_database_file,
                                 SQLITE_OPEN_READWRITE);

      LatticeParameters parameters{.latconst = 1,
                              .boxxhi = 100,
                              .boxyhi = 100,
                              .boxzhi = 2,
                              .temperature = 300,
                              .g_e = -2.1, .is_add_sites = true,
                              .charge_transfer_style = ChargeTransferStyle::MARCUS,
                              .lattice_fill = "none"}; 


      dynamic_LGMC_ = LatticeReactionNetwork(model_database,
                                          initial_state_database,
                                          parameters);

   }

   LatticeReactionNetwork static_LGMC_
   LatticeReactionNetwork dynamic_LGMC_
}

TEST_F(LatticeReactionNetworkTest, InitializeMembersStatic) {

   

}

TEST_F(LatticeReactionNetworkTest, InitializeMembersDynamic) {

   

}




