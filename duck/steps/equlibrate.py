import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import pickle
from duck.utils import duck_stuff
from duck.utils import cal_ints

def do_equlibrate(force_constant_equilibrate=1.0,gpu_id=0):
    # Find the interations
    keyInteraction = cal_ints.find_interaction()
    # Platform definition
    platform = mm.Platform_getPlatformByName("OpenCL")
    platformProperties = {}
    platformProperties['OpenCLPrecision'] = 'mixed'
    platformProperties["OpenCLDeviceIndex"] = gpu_id
    print("loading pickle")
    pickle_in=open('complex_system.pickle', 'rb')
    combined_pmd = pickle.load(pickle_in)[0]
    combined_pmd.symmetry=None
    pickle_in.close()
    ##################
    ##################
    #  Minimisation  #
    ##################
    ##################
    print('Minimising...')
    # Define system
    system = combined_pmd.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom)
    # Get indexes of heavy atoms in chunk
    Chunk_Heavy_Atoms = duck_stuff.getHeavyAtomsInSystem(combined_pmd)
    # Apply force on all havy atoms of chunk
    duck_stuff.applyHarmonicPositionalRestraints(system, force_constant_equilibrate, combined_pmd.positions, Chunk_Heavy_Atoms)
    # Integrator
    integrator = mm.VerletIntegrator(1*u.femtosecond)
    # Define Simulation
    simulation = app.Simulation(combined_pmd.topology, system, integrator, platform,platformProperties)
    simulation.context.setPositions(combined_pmd.positions)
    print(simulation.context.getPlatform().getName())
    for key in simulation.context.getPlatform().getPropertyNames():
        print(key, simulation.context.getPlatform().getPropertyValue(simulation.context, key))
    # Minimizing energy
    simulation.minimizeEnergy(maxIterations=1000)
    # Saving minimised positions
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))
    ##########################
    ##########################
    # Equlibration - heating #
    ##########################
    ##########################
    #new minimised positions, however using old restraints
    # Define new system
    system = combined_pmd.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom, constraints=app.HBonds, hydrogenMass=None)
    # Apply force on all havy atoms of chunk and apply restraint for the ligand-chunk distance
    duck_stuff.applyHarmonicPositionalRestraints(system, force_constant_equilibrate, combined_pmd.positions, Chunk_Heavy_Atoms)
    duck_stuff.applyLigandChunkRestraint(system, force_constant_equilibrate, 10.0, 2*u.angstrom, 3*u.angstrom, 4*u.angstrom, keyInteraction)
    # Intergator
    integrator = mm.LangevinIntegrator(300*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)
    # Define Simulation
    simulation = app.Simulation(combined_pmd.topology, system, integrator, platform,platformProperties)
    simulation.context.setPositions(positions) #changing coordintes to minimized
    # Reporters
    simulation.reporters.append(app.StateDataReporter("heating.csv", 1000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=50000))
    simulation.reporters.append(app.DCDReporter("heating.dcd", 1000))
    # Heating the system
    print("Heating ... ")
    simulation.step(50000) # 0.01 ns
    # Save the positions and velocities
    positions = simulation.context.getState(getPositions=True).getPositions()
    velocities = simulation.context.getState(getVelocities=True).getVelocities()
    app.PDBFile.writeFile(simulation.topology, positions, open('heating_final.pdb', 'w'))
    #clear reporters
    simulation.reporters = []
    ##########################
    ##########################
    # Equlibration - density #
    ##########################
    ##########################
    simulation = duck_stuff.setUpNPTEquilibration(system, combined_pmd,platform, platformProperties, positions, velocities)
    # Reporters
    simulation.reporters.append(app.StateDataReporter("density.csv", 1000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=50000))
    # Correcting the density
    print("Correcting density")
    simulation.step(50000) # 0.01 ns
    # Save the positions and velocities
    positions = simulation.context.getState(getPositions=True).getPositions()
    velocities = simulation.context.getState(getVelocities=True).getVelocities()
    app.PDBFile.writeFile(simulation.topology, positions, open('density_final.pdb', 'w'))
    #saving simulation stage
    checkpoint = 'equil.chk'
    simulation.saveCheckpoint(checkpoint)
    return [checkpoint]

if __name__ == "__main__":
    do_equlibrate()