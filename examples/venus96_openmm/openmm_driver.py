import argparse

import numpy as np
from openmm import *
from openmm.app import *
from openmm.unit import *


def _parse_args():
    parser = argparse.ArgumentParser(
        prog="openmm_driver.py",
        description="A VENUS96 OpenMM Driver.",
    )

    parser.add_argument(
        "-t",
        type=str,
        help="Type: energy, forces",
        required=True,
    )

    parser.add_argument(
        "-f",
        type=str,
        help="Input file: positions.txt",
        required=True,
    )

    args = parser.parse_args()

    return args


def _build_simulation() -> Simulation:
    gro = GromacsGroFile("input.gro")

    top = GromacsTopFile(
        "topol.top",
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
    )

    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
    )

    integrator = LangevinMiddleIntegrator(
        290.0 * kelvin,
        1.0 / picosecond,
        1.0e-03 * picosecond,
    )

    simulation = Simulation(
        topology=top.topology,
        system=system,
        integrator=integrator,
    )

    return simulation


def _load_positions(file: str) -> np.array:
    return np.loadtxt(file).reshape((-1, 3))


if __name__ == "__main__":
    # parse args
    args = _parse_args()

    # build simulation
    simulation = _build_simulation()

    # load and set positions
    positions = _load_positions(args.f) * 0.1  # Angstrom -> nm
    simulation.context.setPositions(positions)

    if args.t == "energy":
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        _energy = energy.value_in_unit(kilocalorie_per_mole)

        with open("energy.txt", "w") as f:
            f.write(f"{_energy:15.9f}\n")
            f.flush()

    elif args.t == "forces":
        state = simulation.context.getState(getForces=True)
        forces = state.getForces()
        _forces = [
            force.value_in_unit(kilocalorie_per_mole / angstrom) for force in forces
        ]

        with open("forces.txt", "w") as f:
            for _force in _forces:
                f.write(
                    f"{_force[0]:15.9f}{_force[1]:15.9f}{_force[2]:15.9f}\n"
                )

            f.flush()
