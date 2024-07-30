/* ----------------------------------------------------------------------
Unit tests for the GMC solvers
All tests use googletest unit test framework
---------------------------------------------------------------------- */

#include <gtest/gtest.h>
#include "../GMC/linear_solver.h"
#include "../GMC/tree_solver.h"
#include "../GMC/sparse_solver.h"

TEST(GMC_solvers, GMC_solvers)
{

    std::vector<double> initial_propensities = {
        0, 0, 0, 0.1,
        0, 0, 0, 0, 0.2,
        0, 0.3,
        0, 0, 0.1,
        0.1, 0, 0};

    TreeSolver tree_solver(42, std::ref(initial_propensities));
    LinearSolver linear_solver_unused(42, std::ref(initial_propensities));
    SparseSolver sparse_solver(42, std::ref(initial_propensities));
    LinearSolver linear_solver(42, std::move(initial_propensities));

    for (int i = 0; i < 100000; i++)
    {
        Event linear_event = linear_solver.event().value();
        Event tree_event = tree_solver.event().value();
        Event sparse_event = sparse_solver.event().value();
        EXPECT_EQ(linear_event.index, tree_event.index);
        EXPECT_EQ(linear_event.index, sparse_event.index);
    }
}