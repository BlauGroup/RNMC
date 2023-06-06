/* -------------------------------------------------------------------------
*
* Unit tests for reaction network
* All tests use catch2 unit test framework
*
------------------------------------------------------------------------- */

#include "../core/sql.h"
#include "../GMC/reaction_network.h"

#include "gtest/gtest.h"
#include <string>

// Using gtest fixture to use same ReactionNetwork 

class ReactionNetworkTest : public ::testing::Test {
   protected:
   void SetUp() override {
      std::string model_database_file = "test_sqlite_files/GMC/network.sqlite";
      std::string initial_state_database_file = "test_sqlite_files/GMC/state.sqlite";

      SqlConnection model_database = SqlConnection(model_database_file,
                                    SQLITE_OPEN_READWRITE);
      SqlConnection initial_state_database = SqlConnection(initial_state_database_file,
                                    SQLITE_OPEN_READWRITE);

      ReactionNetworkParameters parameters;

      reaction_network_ = ReactionNetwork(model_database,
                                          initial_state_database,
                                          parameters);
   }

   ReactionNetwork reaction_network_;

};

TEST_F(ReactionNetworkTest, InitializeMembers) {
   
   // check basic member variables
   EXPECT_EQ(reaction_network_.reactions.size(), 7);
   EXPECT_EQ(reaction_network_.initial_state.size(), 7);
   EXPECT_EQ(reaction_network_.initial_propensities.size(), 7);
   EXPECT_EQ(reaction_network_.factor_zero, 1);
   EXPECT_EQ(reaction_network_.factor_two, 1);
   EXPECT_EQ(reaction_network_.factor_duplicate, 0.5);

}

TEST_F(ReactionNetworkTest, InitializeReactions) {
   
   // one reactant, one product
   EXPECT_EQ(reaction_network_.reactions[0].number_of_reactants, 1);
   EXPECT_EQ(reaction_network_.reactions[0].number_of_products, 1);
   EXPECT_EQ(reaction_network_.reactions[0].reactants[0], 0);
   EXPECT_EQ(reaction_network_.reactions[0].products[0], 2);
   EXPECT_EQ(reaction_network_.reactions[0].rate, 2);

   // two reactants, one product
   EXPECT_EQ(reaction_network_.reactions[1].number_of_reactants, 2);
   EXPECT_EQ(reaction_network_.reactions[1].number_of_products, 1);
   EXPECT_EQ(reaction_network_.reactions[1].reactants[0], 0);
   EXPECT_EQ(reaction_network_.reactions[1].reactants[1], 1);
   EXPECT_EQ(reaction_network_.reactions[1].products[0], 3);
   EXPECT_EQ(reaction_network_.reactions[1].rate, 4);

   // one reactant, two products
   EXPECT_EQ(reaction_network_.reactions[2].number_of_reactants, 1);
   EXPECT_EQ(reaction_network_.reactions[2].number_of_products, 2);
   EXPECT_EQ(reaction_network_.reactions[2].reactants[0], 1);
   EXPECT_EQ(reaction_network_.reactions[2].products[0], 2);
   EXPECT_EQ(reaction_network_.reactions[2].products[1], 3);
   EXPECT_EQ(reaction_network_.reactions[2].rate, 5);

   // two reactants, two products
   EXPECT_EQ(reaction_network_.reactions[3].number_of_reactants, 2);
   EXPECT_EQ(reaction_network_.reactions[3].number_of_products, 2);
   EXPECT_EQ(reaction_network_.reactions[3].reactants[0], 2);
   EXPECT_EQ(reaction_network_.reactions[3].reactants[1], 1);
   EXPECT_EQ(reaction_network_.reactions[3].products[0], 3);
   EXPECT_EQ(reaction_network_.reactions[3].products[1], 4);
   EXPECT_EQ(reaction_network_.reactions[3].rate, 3);
}


TEST_F(ReactionNetworkTest, InitializeState) {
   
   EXPECT_EQ(reaction_network_.initial_state[0], 0);
   EXPECT_EQ(reaction_network_.initial_state[1], 10000);
   EXPECT_EQ(reaction_network_.initial_state[2], 200);
   EXPECT_EQ(reaction_network_.initial_state[3], 20);
   EXPECT_EQ(reaction_network_.initial_state[4], 500);
   EXPECT_EQ(reaction_network_.initial_state[5], 200);
   EXPECT_EQ(reaction_network_.initial_state[6], 100);

}

TEST_F(ReactionNetworkTest, ComputeDependents) {
   EXPECT_EQ(reaction_network_.dependents.size(), 7);
   
   EXPECT_EQ(reaction_network_.dependents[0].size(), 3);
   EXPECT_EQ(reaction_network_.dependents[0][0], 0);
   EXPECT_EQ(reaction_network_.dependents[0][1], 1);
   EXPECT_EQ(reaction_network_.dependents[0][2], 5);


   EXPECT_EQ(reaction_network_.dependents[1].size(), 3);
   EXPECT_EQ(reaction_network_.dependents[1][0], 1);
   EXPECT_EQ(reaction_network_.dependents[1][1], 2);
   EXPECT_EQ(reaction_network_.dependents[1][2], 3);
   
   EXPECT_EQ(reaction_network_.dependents[2].size(), 1);
   EXPECT_EQ(reaction_network_.dependents[2][0], 3);
   
   EXPECT_EQ(reaction_network_.dependents[3].size(), 2);
   EXPECT_EQ(reaction_network_.dependents[3][0], 5);
   EXPECT_EQ(reaction_network_.dependents[3][1], 6);

   EXPECT_EQ(reaction_network_.dependents[4].size(), 0);
   EXPECT_EQ(reaction_network_.dependents[5].size(), 1);
   EXPECT_EQ(reaction_network_.dependents[5][0], 4);

   EXPECT_EQ(reaction_network_.dependents[6].size(), 1);
   EXPECT_EQ(reaction_network_.dependents[6][0], 4);

}


TEST_F(ReactionNetworkTest, InitializePropensities) {
   
   for(int i = 0; i < 7; i++) {
      EXPECT_EQ(reaction_network_.initial_propensities[i], 0);
   }
   
}

// compute propensity 
// update_state 
// update_propensities
// checkpoint
// store_checkpoint

