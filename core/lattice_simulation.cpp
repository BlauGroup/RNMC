#include "lattice_simulation.h"

void LatticeSimulation::init() {

    lattice_network.compute_initial_propensities(state.homogeneous, state.lattice);

    latSolver = LatticeSolver(seed, std::ref(lattice_network.initial_propensities));
    this->update_function = [&] (Update update) {latSolver.update(update);};
    lattice_update_function = [&] (LatticeUpdate lattice_update, 
                                    std::unordered_map<std::string,
                                    std::vector< std::pair<double, int> > > &props) 
                                    {latSolver.update(lattice_update, props);};

    lattice_network.update_adsorp_state(state.lattice, this->props, 
                                        latSolver.propensity_sum, 
                                        latSolver.number_of_active_indices);
    lattice_network.update_adsorp_props(state.lattice, lattice_update_function, 
                                        state.homogeneous, std::ref(props));

    // only call if checkpointing 
    if(state.lattice->isCheckpoint) {
        lattice_network.update_all_propensities(state.lattice, props, 
                                                latSolver.propensity_sum, 
                                                latSolver.number_of_active_indices, 
                                                lattice_update_function);
        lattice_network.print_state_propensities(latSolver.propensity_sum, 
                                                 latSolver.propensities, 
                                                 state.homogeneous, "checkpoint.txt");     
    }
} //init()

/*---------------------------------------------------------------------------*/

bool LatticeSimulation::execute_step() {

    std::optional<LatticeEvent> maybe_event = latSolver.event_lattice(props);

    if (!maybe_event) {

        return false;

    } else {
        // an event happens
        LatticeEvent event = maybe_event.value();
        int next_reaction = event.index;

        // randomly assign sites when two products
        // use boolean to determine order for trajectories
        bool flip_sites = false;

        // update_state
        lattice_network.update_state(state.lattice, std::ref(props), 
                    std::ref(this->state.homogeneous), next_reaction, 
                    event.site_one, event.site_two, latSolver.propensity_sum, 
                    latSolver.number_of_active_indices, flip_sites);

        if(state.lattice->sites.size() != state.lattice->idneigh.size()) {
            assert(false);
        }
        if(state.lattice->loc_map.size() != state.lattice->idneigh.size()) {
            assert(false);
        }
        // update time
        this->time += event.dt;
        int site_1_mapping;
        int site_2_mapping;

        if(event.site_one) {
            int site_1 = event.site_one.value();

            if(site_1 == -2) {
                site_1_mapping = -2;
            }
            else if (site_1 == -3){
                site_1_mapping = -3;
            }
            else {
                assert(state.lattice->sites.find(site_1) != state.lattice->sites.end());
                site_1_mapping = lattice_network.combine(state.lattice->sites[site_1].i, 
                                                         state.lattice->sites[site_1].j, 
                                                         state.lattice->sites[site_1].k);
            }
        }
        if(event.site_two) {
            int site_2 = event.site_two.value();

            if(site_2 == -2) {
                site_2_mapping = -2;
            }
            else if (site_2 == -3){
                site_2_mapping = -3;
            }
            else {
                assert(state.lattice->sites.find(site_2) != state.lattice->sites.end());
                site_2_mapping = lattice_network.combine(state.lattice->sites[site_2].i, 
                                                         state.lattice->sites[site_2].j, 
                                                         state.lattice->sites[site_2].k);
            }
        }
        if(flip_sites) {
            // order of sites must corrrespond to order of products in reaction
            std::swap(site_1_mapping, site_2_mapping);
        }

        // record what happened
        history.push_back(LatticeTrajectoryHistoryElement {
            .seed = this->seed,
            .reaction_id = next_reaction,
            .time = this->time,
            .step = this->step,
            .site_1_mapping = site_1_mapping,
            .site_2_mapping = site_2_mapping
            });

        if (history.size() == this->history_chunk_size ) {
            history_queue.insert_history(
                HistoryPacket<LatticeTrajectoryHistoryElement> {
                    .history = std::move(this->history),
                    .seed = this->seed
                    });

            history = std::vector<LatticeTrajectoryHistoryElement> ();
            history.reserve(this->history_chunk_size);
        }

        // update_propensities 
        lattice_network.update_propensities(state.lattice, 
                        std::ref(this->state.homogeneous), 
                        this->update_function, lattice_update_function, 
                        next_reaction, event.site_one, event.site_two, props);

        // increment step
        this->step++;

        return true;
    }
 
} //execute_step()


/*---------------------------------------------------------------------------*/

void LatticeSimulation::print_output() {
    lattice_network.print_state_propensities(latSolver.propensity_sum, 
                                             latSolver.propensities, 
                                             state.homogeneous, "final_state.txt");
} // print_output()