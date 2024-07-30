/* ----------------------------------------------------------------------
Unit tests for nano particle
All tests use googletest unit test framework
---------------------------------------------------------------------- */

#include "gtest/gtest.h"
#include <string>

#include "../core/sql.h"
#include "../NPMC/nano_particle.h"

class NanoParticleTEST : public ::testing::Test
{
protected:
    void SetUp() override
    {

        std::string model_database_file = "../examples/NPMC/simple_example/np.sqlite";
        std::string initial_state_database_file = "../examples/NPMC/simple_example/initial_state.sqlite";

        SqlConnection model_database = SqlConnection(model_database_file,
                                                     SQLITE_OPEN_READWRITE);
        SqlConnection initial_state_database = SqlConnection(initial_state_database_file,
                                                             SQLITE_OPEN_READWRITE);

        NanoParticleParameters parameters;

        nano_particle_ = NanoParticle(model_database,
                                      initial_state_database,
                                      parameters);
    }

    NanoParticle nano_particle_;
};

TEST_F(NanoParticleTEST, InitializeDegreesOfFreedom)
{
    EXPECT_EQ(nano_particle_.degrees_of_freedom[0], 2);
    EXPECT_EQ(nano_particle_.degrees_of_freedom[1], 4);
    EXPECT_EQ(nano_particle_.degrees_of_freedom[2], 5);
}

TEST_F(NanoParticleTEST, InitializeSites)
{
    EXPECT_EQ(static_cast<int>(nano_particle_.sites.size()), 5);

    std::vector<double> x_cords{-0.30331966499999991, 0., 0., -0.30335000000000012, 0.};
    std::vector<double> y_cords{0.17512169023825716, 0., 0., -0.52541761247601904, 0.};
    std::vector<double> z_cords{0.17757499999999987, -0.35294807000000006, 0.35735192999999993, -0.0022019300000001075, -0.0022019300000000186};
    std::vector<double> species_id{0, 1, 1, 2, 2};

    for (int i = 0; i < static_cast<int>(nano_particle_.sites.size()); i++)
    {
        EXPECT_EQ(nano_particle_.sites[i].x, x_cords[i]);
        EXPECT_EQ(nano_particle_.sites[i].y, y_cords[i]);
        EXPECT_EQ(nano_particle_.sites[i].z, z_cords[i]);
        EXPECT_EQ(nano_particle_.sites[i].species_id, species_id[i]);
    }
}

TEST_F(NanoParticleTEST, SiteDistanceSquared)
{

    NanoSite s1 = {0, 0, 0, 0};
    NanoSite s2 = {0, 1, 1, 1};
    NanoSite s3 = {0, 0, 0, 2};
    NanoSite s4 = {0.5, 0.3, 1.5, 3};

    EXPECT_EQ(nano_particle_.site_distance_squared(s1, s2), 2);
    EXPECT_EQ(nano_particle_.site_distance_squared(s1, s3), 0);
    EXPECT_EQ(nano_particle_.site_distance_squared(s2, s4), 0.99);
}

TEST_F(NanoParticleTEST, DistanceMatrix)
{
    // SiteDistanceSquare previously tested
    // testing location of distances in distance matrix
    for (int i = 0; i < static_cast<int>(nano_particle_.distance_matrix.size()); i++)
    {
        for (int j = 0; j < static_cast<int>(nano_particle_.distance_matrix[i].size()); j++)
        {
            EXPECT_EQ(nano_particle_.distance_matrix[i][j],
                      std::sqrt(nano_particle_.site_distance_squared(nano_particle_.sites[i], nano_particle_.sites[j])));
        }
    }
}
