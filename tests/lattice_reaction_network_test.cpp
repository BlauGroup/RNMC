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

// check initalization
TEST(create_lattice_reaction_network, LatticeReactionNetwork) {

   std::cout << "I am here!" << std::endl;

   std::string model_database_file = "../examples/LGMC/CO_oxidation(static)/network.sqlite";
   std::string initial_state_database_file = "../examples/LGMC/CO_oxidation(static)/state.sqlite";

   std::cout << "I am here!" << std::endl;

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


   LatticeReactionNetwork model = LatticeReactionNetwork(model_database,
                                       initial_state_database,
                                       parameters);


}


// correctly initalizes dependents
// correctly computes propensities
// check


