#include "../core/sql.h"
#include "../NPMC/nano_particle.h"

#include "gtest/gtest.h"
#include <string>

// Using gtest fixture to use same NanoParticle 

class NanoParticleTEST : public ::testing::Test {
   protected:
   void SetUp() override {

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

TEST_F(NanoParticleTEST, InitializeMembers) {
    EXPECT_EQ(static_cast<int> (nano_particle_.sites.size()), 5);
}