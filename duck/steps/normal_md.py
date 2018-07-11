import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import sys, os
import pickle
from duck.utils import duck_stuff
from duck.utils import cal_ints


def perform_md(checkpoint_in_file, checkpoint_out_file, csv_out_file, pdb_out_file, force_constant_ligand=1.0,
               md_len=1.0, force_constant_chunk=0.1, gpu_id=0):
    if os.path.isfile(checkpoint_out_file):
        return
    print("loading pickle")
    pickle_in = open('complex_system.pickle', 'rb')
    combined_pmd = pickle.load(pickle_in)[0]
    print(dir(combined_pmd))
    key_interaction = cal_ints.find_interaction()
    pickle_in.close()
    traj_out_file = "traj.out"
    MD_len = md_len * u.nanosecond
    sim_steps = round(MD_len / (0.002 * u.picosecond))
    # Platform definition
    platform = mm.Platform_getPlatformByName("OpenCL")
    platformProperties = {}
    platformProperties['OpenCLPrecision'] = 'mixed'
    platformProperties["OpenCLDeviceIndex"] = gpu_id
    # Get indexes of heavy atoms in chunk
    Chunk_Heavy_Atoms = duck_stuff.getHeavyAtomsInSystem(combined_pmd)
    # Setting System
    system = combined_pmd.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9 * u.angstrom, constraints=app.HBonds,
                                       hydrogenMass=None)
    # Apply force on all havy atoms of chunk and apply restraint for the ligand-chunk distance
    duck_stuff.applyHarmonicPositionalRestraints(system, force_constant_chunk, combined_pmd.positions,
                                                 Chunk_Heavy_Atoms)
    duck_stuff.applyLigandChunkRestraint(system, force_constant_ligand, 10.0, 2 * u.angstrom, 3 * u.angstrom,
                                         4 * u.angstrom, key_interaction)
    # Integrator
    integrator = mm.LangevinIntegrator(300 * u.kelvin, 4 / u.picosecond, 0.002 * u.picosecond)
    # Setting Simulation object and loading the checkpoint
    simulation = app.Simulation(combined_pmd.topology, system, integrator, platform, platformProperties)
    simulation.loadCheckpoint(checkpoint_in_file)
    # Simulation reporters
    simulation.reporters.append(
        app.StateDataReporter(csv_out_file, 2000, step=True, time=True, totalEnergy=True, kineticEnergy=True,
                              potentialEnergy=True, temperature=True, density=True, progress=True, totalSteps=sim_steps,
                              speed=True))
    simulation.reporters.append(app.DCDReporter("md.dcd", 100000))
    # Production
    simulation.step(sim_steps)
    # Save state in checkpoint file and save coordinates in PDB file
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    positions = state.getPositions()
    velocities = state.getVelocities()
    app.PDBFile.writeFile(simulation.topology, positions, open(pdb_out_file, 'w'))
    simulation.saveCheckpoint(checkpoint_out_file)


if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit("Usage 03_md.py in.chk out.chk out.csv out.pdb")
    checkpoint_in_file = sys.argv[1]
    checkpoint_out_file = sys.argv[2]
    csv_out_file = sys.argv[3]
    pdb_out_file = sys.argv[4]
    perform_md(checkpoint_in_file, checkpoint_out_file, csv_out_file, pdb_out_file)
