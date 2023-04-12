import os
import sys
import subprocess
import pickle

from HiPRGen.network_loader import NetworkLoader
from HiPRGen.initial_state import find_mol_entry_from_xyz_and_charge
from HiPRGen.report_generator import ReportGenerator
from HiPRGen.initial_state import insert_initial_state

from rnmcpy.analysis.gmc_analysis import (
    reaction_tally_report,
    species_report,
    Pathfinding,
    SimulationReplayer,
    generate_pathway_report,
    sink_report,
    consumption_report,
)

# Since HiPRGen uses an end-to-end testing approach rather than testing
# each individual function, we have decided to use the tests as
# documentation, by explaining every single line through the first test.


# The first thing you need to consider when using HiPRGen is how many
# worker threads do you want to run. HiPRGen can be run with a single
# thread or thousands distrubuted across several nodes. For reaction
# networks with between ~5000 and ~10000 species, we have found that the
# optimal number of worker threads is between 1000 and 2000. If you try
# and use more than that, the worker threads are going to spend lots of
# time waiting for the dispatcher to get through all of the reactions it
# is being sent, which slows everything down. Fixing this would require
# a more complex distrubuted system, but it hasn't been an issue even
# for very large reaction networks.
if len(sys.argv) != 2:
    print("usage: python test.py number_of_threads")
    quit()


number_of_threads = sys.argv[1]


class bcolors:
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'


# HiPRGen is organized as a pipeline, where all the relevent data is
# stored in a sqlite database between phases. For this reason, during
# a single run of the full pipeline, it makes sense to store all the
# relevent files in a single directory. We have two test sets, a lithium
# set and a magnesium set. Since the lithium test set is older, we shall
# document that instead of the mg test set.

if os.path.isdir('./scratch'):
    subprocess.run(['rm', '-r', './scratch'])

subprocess.run(['mkdir', './scratch'])


def li_test():
    # folder is the where we store all our intermediate databases
    folder = './scratch/li_test'
    subprocess.run(['mkdir', folder])

    # There is one non-local part of species filtering: we consider two
    # molecules to be equivalent if they have the same total charge,
    # composition, and covalent bonds, even if they have different metal
    # coordination, and we choose one such molecule in each "coordimer"
    # class using the coodimer weight function.  Since most of our logging
    # later on is defined in terms of a fixed molecule set, logging for
    # the species filtering phase is messy, so ignore the species_report
    # argument for now. The second argument is where we store a pickle of
    # the filtered molecule entries for use in later phases.

    mol_entries = pickle.load("../data/li_mol_entries.pickle")

    # After we have generated the mol_entries, we refer to molecules by
    # their index. The function find_mol_entry_from_xyz_and_charge can
    # help find the indices of specific species to be used in the initial
    # condition for propagating trajectories and/or trajectory analysis.
    Li_plus_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        '../test_materials/xyz_files/Li.xyz',
        1)

    EC_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        '../test_materials/xyz_files/EC.xyz',
        0)

    LEDC_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        '../test_materials/xyz_files/LEDC.xyz',
        0)

    # After generating a reaction network, it is stored in rn.sqlite. We
    # use Monte Carlo simulation to interrogate the network, and for that
    # we need to define an initial condition.
    initial_state = {
        Li_plus_id: 30,
        EC_id: 30
    }

    # The initial state and the trajectories (after simulation) are stored in
    # a seperate database from the network, here called initial_state.sqlite.
    # This facilitates running multiple independent simulations of the same
    # network with different initial conditions at the same time, if desired.
    insert_initial_state(initial_state, mol_entries, folder + '/initial_state.sqlite')

    subprocess.run([
        'GMC',
        '--reaction_database=' + '../data/li_rn.sqlite',
        '--initial_state_database=' + folder + '/initial_state.sqlite',
        '--number_of_simulations=1000',
        '--base_seed=1000',
        '--thread_count=' + number_of_threads,
        '--step_cutoff=200'
    ])

    # The network loader builds a python object around a reaction network
    # and the molecules to make it easier to use them.
    network_loader = NetworkLoader(
        '../data/li_rn.sqlite',
        '../data/li_mol_entries.pickle',
        folder + '/initial_state.sqlite'
        )

    network_loader.load_trajectories()
    network_loader.load_initial_state()

    # HiPRGen has analysis tools to understand what happened in our simulation.
    # The output files are written into the same folder in which the reaction
    # network is stored.

    # This report is empty, but we use it to generate the molecule pictures.
    # This is an expensive operation, so we only want do do it once.
    _ = ReportGenerator(
        network_loader.mol_entries,
        folder + '/dummy.tex',
        rebuild_mol_pictures=True)

    # The tally report shows reactions sorted by the number of times fired.
    reaction_tally_report(
        network_loader,
        folder + '/reaction_tally.tex'
    )
    # Run `pdflatex reaction_tally.tex` in `scratch/li_test` to generate
    # the tally report PDF.

    # The species report shows every specie in the network and their IDs.
    species_report(network_loader, folder + '/species_report.tex')
    # Run `pdflatex species_report.tex` in `scratch/li_test` to generate
    # the species report PDF.

    # Pathfinding is a central goal of HiPRGen / GMC. See mc_analysis.py for
    # further documentation of the Pathfinding class.
    pathfinding = Pathfinding(network_loader)

    # The pathway report shows all the ways that a target species was
    # produced in the simulation trajectories, where each simulation only
    # contributes the shortest path responsible for the first formation
    # of the target species to the report. The report can be sorted by
    # pathway frequency, but instead here we sort by pathway cost. Note
    # that the test network has ~5000 reactions while production networks
    # have between 50-100 million reactions.
    generate_pathway_report(
        pathfinding,
        LEDC_id,
        folder + '/LEDC_pathways.tex',
        sort_by_frequency=False
    )
    # Run `pdflatex LEDC_pathways.tex` in `scratch/li_test` to generate
    # the LEDC pathway report PDF.

    # The simulation replayer sweeps through all trajectories in order
    # to extract additional information that is used for consumption
    # reports and sink reports.
    simulation_replayer = SimulationReplayer(network_loader)

    # The consumption report shows reactions which consumed a target
    # species, sorted by the number of times the reaction fired.
    consumption_report(simulation_replayer,
                       LEDC_id,
                       folder + '/LEDC_consumption_report.tex')
    # Run `pdflatex LEDC_consumption_report.tex` in `scratch/li_test`
    # to generate the LEDC consumption report PDF.

    # The sink report shows species which have a production to
    # consumption ratio of greater than 3/2 and which have an expected
    # value above 0.1. These are two of the three heuristic criteria
    # that we use to identify network products. The third criteria is
    # that each network product must have a shortest path with cost
    # less than 10. This can be checked by generating pathway reports
    # to each species shown in the sink report. For the curious reader,
    # we note that generating pathway reports to the six species in the
    # sink report will show that only Li2CO3, C2H4, LiEDC-, and DLEMC
    # have sufficiently low-cost paths to pass the third criteria and
    # thus to be considered products of the test network used here.
    sink_report(simulation_replayer, folder + '/sink_report.tex')
    # Run `pdflatex sink_report.tex` in `scratch/li_test` to generate
    # the sink report PDF.

    return True


def mg_test():
    folder = './scratch/mg_test'
    subprocess.run(['mkdir', folder])

    mol_entries = pickle.load("../data/mg_mol_entries.pickle")

    c2h4_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        '../test_materials/xyz_files/c2h4.xyz',
        0)

    c2h6_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        '../test_materials/xyz_files/c2h6.xyz',
        0)

    initial_state = {
        33: 30,
        81: 30
    }

    insert_initial_state(initial_state, mol_entries, folder + '/initial_state.sqlite')

    subprocess.run([
        'GMC',
        '--reaction_database=' + '../data/mg_rn.sqlite',
        '--initial_state_database=' + folder + '/initial_state.sqlite',
        '--number_of_simulations=1000',
        '--base_seed=1000',
        '--thread_count=' + number_of_threads,
        '--step_cutoff=200'
    ])

    network_loader = NetworkLoader(
        '../data/mg_rn.sqlite',
        '../data/mg_mol_entries.pickle',
        folder + '/initial_state.sqlite'
        )

    network_loader.load_trajectories()
    network_loader.load_initial_state()

    _ = ReportGenerator(
        network_loader.mol_entries,
        folder + '/dummy.tex',
        rebuild_mol_pictures=True)

    reaction_tally_report(
        network_loader,
        folder + '/reaction_tally.tex'
    )

    pathfinding = Pathfinding(network_loader)

    generate_pathway_report(
        pathfinding,
        c2h6_id,
        folder + '/C2H6_pathways.tex',
        sort_by_frequency=False
    )

    generate_pathway_report(
        pathfinding,
        c2h4_id,
        folder + '/C2H4_pathways.tex',
        sort_by_frequency=False
    )

    species_report(network_loader, folder + '/species_report.tex')

    return True


tests = [
    mg_test,
    li_test,
]

for test in tests:
    if not test():
        exit(1)
