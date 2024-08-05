/* ----------------------------------------------------------------------
Unit tests for reaction network
All tests use googletest unit test framework
---------------------------------------------------------------------- */

#include <string>

#include "../core/sql.h"
#include "../GMC/gillespie_reaction_network.h"
#include "../GMC/tree_solver.h"
#include "gtest/gtest.h"

class ReactionNetworkTest : public ::testing::Test
{
protected:
   void SetUp() override
   {
      std::string model_database_file = "test_sqlite_files/GMC/network.sqlite";
      std::string initial_state_database_file = "test_sqlite_files/GMC/state.sqlite";

      SqlConnection model_database = SqlConnection(model_database_file,
                                                   SQLITE_OPEN_READWRITE);
      SqlConnection initial_state_database = SqlConnection(initial_state_database_file,
                                                           SQLITE_OPEN_READWRITE);

      ReactionNetworkParameters parameters;

      reaction_network_ = GillespieReactionNetwork(model_database,
                                                   initial_state_database,
                                                   parameters);
      
      reaction_network_.compute_initial_propensities(reaction_network_.initial_state, initial_props);
   }

   GillespieReactionNetwork reaction_network_;
   std::vector<double> initial_props;
};

TEST_F(ReactionNetworkTest, InitializeMembers)
{

   // check basic member variables
   EXPECT_EQ(int(reaction_network_.reactions.size()), 7);
   EXPECT_EQ(int(reaction_network_.initial_state.size()), 7);
   EXPECT_EQ(reaction_network_.factor_zero, 1);
   EXPECT_EQ(reaction_network_.factor_two, 1);
   EXPECT_EQ(reaction_network_.factor_duplicate, 0.5);
}

TEST_F(ReactionNetworkTest, InitializeReactions)
{

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

TEST_F(ReactionNetworkTest, InitializeState)
{

   EXPECT_EQ(reaction_network_.initial_state[0], 0);
   EXPECT_EQ(reaction_network_.initial_state[1], 10000);
   EXPECT_EQ(reaction_network_.initial_state[2], 200);
   EXPECT_EQ(reaction_network_.initial_state[3], 20);
   EXPECT_EQ(reaction_network_.initial_state[4], 500);
   EXPECT_EQ(reaction_network_.initial_state[5], 200);
   EXPECT_EQ(reaction_network_.initial_state[6], 100);
}

TEST_F(ReactionNetworkTest, ComputeDependents)
{
   EXPECT_EQ(int(reaction_network_.dependents.size()), 7);

   EXPECT_EQ(int(reaction_network_.dependents[0].size()), 3);
   EXPECT_EQ(reaction_network_.dependents[0][0], 0);
   EXPECT_EQ(reaction_network_.dependents[0][1], 1);
   EXPECT_EQ(reaction_network_.dependents[0][2], 5);

   EXPECT_EQ(int(reaction_network_.dependents[1].size()), 3);
   EXPECT_EQ(reaction_network_.dependents[1][0], 1);
   EXPECT_EQ(reaction_network_.dependents[1][1], 2);
   EXPECT_EQ(reaction_network_.dependents[1][2], 3);

   EXPECT_EQ(int(reaction_network_.dependents[2].size()), 1);
   EXPECT_EQ(reaction_network_.dependents[2][0], 3);

   EXPECT_EQ(int(reaction_network_.dependents[3].size()), 2);
   EXPECT_EQ(reaction_network_.dependents[3][0], 5);
   EXPECT_EQ(reaction_network_.dependents[3][1], 6);

   EXPECT_EQ(int(reaction_network_.dependents[4].size()), 0);

   EXPECT_EQ(int(reaction_network_.dependents[5].size()), 1);
   EXPECT_EQ(reaction_network_.dependents[5][0], 4);

   EXPECT_EQ(int(reaction_network_.dependents[6].size()), 1);
   EXPECT_EQ(reaction_network_.dependents[6][0], 4);
}

TEST_F(ReactionNetworkTest, InitializePropensities)
{

   std::vector<double> expected_propensities = {0, 0, 50000, 6e6, 6e4, 0, 760};
   for (int i = 0; i < 7; i++)
   {
      EXPECT_EQ(initial_props[i], expected_propensities[i]);
   }
}

TEST_F(ReactionNetworkTest, ComputePropensity)
{

   std::vector<double> expected_propensities = {0, 0, 50000, 6e6, 6e4, 0, 760};
   for (int i = 0; i < 7; i++)
   {
      EXPECT_EQ(reaction_network_.compute_propensity(
                    reaction_network_.initial_state, i),
                expected_propensities[i]);
   }
}

TEST_F(ReactionNetworkTest, UpdateState)
{
   std::vector<int> state = {100, 200, 300, 400, 500, 600, 700};
   reaction_network_.update_state(std::ref(state), 0);

   EXPECT_EQ(state[0], 99);
   EXPECT_EQ(state[1], 200);
   EXPECT_EQ(state[2], 301);
   EXPECT_EQ(state[3], 400);
   EXPECT_EQ(state[4], 500);
   EXPECT_EQ(state[5], 600);
   EXPECT_EQ(state[6], 700);

   reaction_network_.update_state(std::ref(state), 1);

   EXPECT_EQ(state[0], 98);
   EXPECT_EQ(state[1], 199);
   EXPECT_EQ(state[2], 301);
   EXPECT_EQ(state[3], 401);
   EXPECT_EQ(state[4], 500);
   EXPECT_EQ(state[5], 600);
   EXPECT_EQ(state[6], 700);

   reaction_network_.update_state(std::ref(state), 2);

   EXPECT_EQ(state[0], 98);
   EXPECT_EQ(state[1], 198);
   EXPECT_EQ(state[2], 302);
   EXPECT_EQ(state[3], 402);
   EXPECT_EQ(state[4], 500);
   EXPECT_EQ(state[5], 600);
   EXPECT_EQ(state[6], 700);

   reaction_network_.update_state(std::ref(state), 3);

   EXPECT_EQ(state[0], 98);
   EXPECT_EQ(state[1], 197);
   EXPECT_EQ(state[2], 301);
   EXPECT_EQ(state[3], 403);
   EXPECT_EQ(state[4], 501);
   EXPECT_EQ(state[5], 600);
   EXPECT_EQ(state[6], 700);

   reaction_network_.update_state(std::ref(state), 4);

   EXPECT_EQ(state[0], 99);
   EXPECT_EQ(state[1], 198);
   EXPECT_EQ(state[2], 301);
   EXPECT_EQ(state[3], 403);
   EXPECT_EQ(state[4], 501);
   EXPECT_EQ(state[5], 599);
   EXPECT_EQ(state[6], 699);

   reaction_network_.update_state(std::ref(state), 5);

   EXPECT_EQ(state[0], 98);
   EXPECT_EQ(state[1], 198);
   EXPECT_EQ(state[2], 302);
   EXPECT_EQ(state[3], 402);
   EXPECT_EQ(state[4], 501);
   EXPECT_EQ(state[5], 599);
   EXPECT_EQ(state[6], 699);

   reaction_network_.update_state(std::ref(state), 6);

   EXPECT_EQ(state[0], 98);
   EXPECT_EQ(state[1], 198);
   EXPECT_EQ(state[2], 303);
   EXPECT_EQ(state[3], 400);
   EXPECT_EQ(state[4], 501);
   EXPECT_EQ(state[5], 599);
   EXPECT_EQ(state[6], 699);
}

TEST_F(ReactionNetworkTest, UpdatePropensities)
{
   TreeSolver tree_solver = TreeSolver(1, initial_props);

   std::vector<int> state = {0, 10000, 200, 20, 500, 200, 100};

   std::function<void(Update)> update_function = [&](Update update)
   { tree_solver.update(update); };
   reaction_network_.update_state(std::ref(state), 4);
   reaction_network_.update_propensities(update_function, std::ref(state), 4);

   EXPECT_EQ(tree_solver.get_propensity(0), 2);
   EXPECT_EQ(tree_solver.get_propensity(1), 40004);
}

// checkpoint
// store_checkpoint
