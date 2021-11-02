Red="\033[0;31m"          # Red
Green="\033[0;32m"        # Green
Color_Off="\033[0m"       # Text Reset

RNMC_TEST_DIR="./test_materials/RNMC"

cp $RNMC_TEST_DIR/initial_state.sqlite $RNMC_TEST_DIR/initial_state_copy.sqlite

./build/RNMC --reaction_database=$RNMC_TEST_DIR/rn.sqlite --initial_state_database=$RNMC_TEST_DIR/initial_state_copy.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200 --dependency_threshold=1

sql='SELECT seed, step, reaction_id FROM trajectories ORDER BY seed ASC, step ASC;'

sqlite3 $RNMC_TEST_DIR/initial_state_with_trajectories.sqlite "${sql}" > $RNMC_TEST_DIR/trajectories
sqlite3 $RNMC_TEST_DIR/initial_state_copy.sqlite "${sql}" > $RNMC_TEST_DIR/copy_trajectories

if  cmp $RNMC_TEST_DIR/trajectories $RNMC_TEST_DIR/copy_trajectories > /dev/null; then
    echo -e "${Green} passed: no difference in trajectories ${Color_Off}"
    RC=0
else
    echo -e "${Red} failed: difference in trajectories ${Color_Off}"
    RC=1
fi

rm $RNMC_TEST_DIR/initial_state_copy.sqlite
rm $RNMC_TEST_DIR/trajectories
rm $RNMC_TEST_DIR/copy_trajectories
exit $RC
