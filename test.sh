Red="\033[0;31m"          # Red
Green="\033[0;32m"        # Green
Color_Off="\033[0m"       # Text Reset

function test_core {
    if ./build/test_core
    then
        echo -e "${Green} passed: linear, sparse and tree samplers agree ${Color_Off}"
        RC=0
    else
        echo -e "${Red} failed: linear, sparse and tree samplers disagree ${Color_Off}"
        RC=1
    fi

}

function test_gmc {
    GMC_TEST_DIR="./test_materials/GMC"

    cp $GMC_TEST_DIR/initial_state.sqlite $GMC_TEST_DIR/initial_state_copy.sqlite

    ./build/GMC --reaction_database=$GMC_TEST_DIR/rn.sqlite --initial_state_database=$GMC_TEST_DIR/initial_state_copy.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=2 --step_cutoff=200 &> /dev/null

    sql='SELECT seed, step, reaction_id FROM trajectories ORDER BY seed ASC, step ASC;'

    sqlite3 $GMC_TEST_DIR/initial_state_with_trajectories.sqlite "${sql}" > $GMC_TEST_DIR/trajectories
    sqlite3 $GMC_TEST_DIR/initial_state_copy.sqlite "${sql}" > $GMC_TEST_DIR/copy_trajectories

    if  cmp $GMC_TEST_DIR/trajectories $GMC_TEST_DIR/copy_trajectories > /dev/null
    then
        echo -e "${Green} passed: no difference in GMC trajectories ${Color_Off}"
        RC=0
    else
        echo -e "${Red} failed: difference in GMC trajectories ${Color_Off}"
        RC=1
    fi

    rm $GMC_TEST_DIR/initial_state_copy.sqlite
    rm $GMC_TEST_DIR/trajectories
    rm $GMC_TEST_DIR/copy_trajectories

}

function test_npmc {
    NPMC_TEST_DIR="./test_materials/NPMC"

    cp $NPMC_TEST_DIR/initial_state.sqlite $NPMC_TEST_DIR/initial_state_copy.sqlite

    # to check for leaks with valgrind, you need to use the option --fair-sched=yes

    ./build/NPMC --nano_particle_database=$NPMC_TEST_DIR/np.sqlite --initial_state_database=$NPMC_TEST_DIR/initial_state_copy.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=2 --step_cutoff=200 # &> /dev/null

    sql='SELECT seed, step, site_id_1, site_id_2, interaction_id FROM trajectories ORDER BY seed ASC, step ASC;'

    sqlite3 $NPMC_TEST_DIR/initial_state_with_trajectories.sqlite "${sql}" > $NPMC_TEST_DIR/trajectories
    sqlite3 $NPMC_TEST_DIR/initial_state_copy.sqlite "${sql}" > $NPMC_TEST_DIR/copy_trajectories

    if  cmp $NPMC_TEST_DIR/trajectories $NPMC_TEST_DIR/copy_trajectories # > /dev/null
    then
        echo -e "${Green} passed: no difference in NPMC trajectories ${Color_Off}"
        RC=0
    else
        echo -e "${Red} failed: difference in NPMC trajectories ${Color_Off}"
        RC=1
    fi

    # rm $NPMC_TEST_DIR/initial_state_copy.sqlite
    # rm $NPMC_TEST_DIR/trajectories
    # rm $NPMC_TEST_DIR/copy_trajectories
}

function check_result {
    if [[ $RC -ne 0 ]]
    then
        exit $RC
    fi
}


test_core
check_result
test_gmc
check_result
test_npmc
check_result

exit $RC
