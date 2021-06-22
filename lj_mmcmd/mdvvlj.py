#!/usr/bin/env python
"""
This module is used for running a simple MD simulations of Lennard-Jones particles
"""

import numpy as np
import pint

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity
from scipy.spatial import distance


class MDvvlj:
    """
    This module is used for running a simple MD simulations of Lennard-Jones particles
    """

    def __init__(self, topology, system_size):
        self.epsilon = Q_(0.238, "kcal/mol").magnitude
        self.sigma = Q_(3.4, ureg.angstrom).magnitude
        self.R = Q_(0.001985875, "kcal/mol/K").magnitude
        self.temperature = Q_(298, ureg.kelvin).magnitude
        self.trajectories = []
        self.potential_energies = []
        self.kinetic_energies = []
        self.deltaT = 2.5e-15  # unit: second
        self.topology = topology
        self.natoms = self.topology.shape[0]
        self.mass = Q_(39.9 * ureg.dalton).to(ureg.kilogram).magnitude
        self.nparticles = self.topology.shape[0]
        self.system_size = system_size
        self.velocities = []

    def _pbcs(self, p1, p2):
        """
        :param p1: position 1

        :param p2: position 2

        :return: distance matrix with periodic boundary conditions
        """
        new = p1 - p2 - self.system_size * np.round((p1 - p2) / self.system_size, 0)
        # return np.sqrt(np.sum(np.square(new), axis=-1))
        return np.linalg.norm(new)

    def pbcs_distance(self, trajectory):
        """
        Calculate distance matrix with periodic boundary conditions

        :param trajectory: coordinates

        :return: this_distance
            distance matrix with periodic boundary conditions
        """

        this_distance = np.zeros((self.nparticles, self.nparticles))
        for atom1 in range(self.nparticles):
            for atom2 in range(self.nparticles):
                this_distance[atom1, atom2] = self._pbcs(trajectory[atom1], trajectory[atom2])
        this_distance = np.tril(this_distance)
        return this_distance

    def calc_potential_energy(self, trajectory):
        """
        Calculate the potential energy of lennard jones fluids

        .. math:: V(r) = 4\epsilon[(\sigma/r)^{12} - (\sigma/r)^6]

        switching function

        .. math:: S = 1 - 6x^5 + 15x^4 - 10x^3

        http://docs.openmm.org/latest/userguide/theory.html#lennard-jones-interaction

        :param trajectory: coordinates

        :return: potential energies matrix
        """
        this_distance = self.pbcs_distance(trajectory)
        shift = lambda x: 1 - 6 * np.power(x, 5) + 15 * np.power(x, 4) - 10 * np.power(x, 3)
        this_distance = np.where(this_distance > 3 * self.sigma, 0.0, this_distance)
        scaled = np.where(this_distance > 2 * self.sigma, shift(((this_distance - 2 * self.sigma) / (self.sigma))), 1.0)
        this_distance = np.where(this_distance == 0.0, 0.0, np.reciprocal(this_distance))

        potential_energy = 4 * self.epsilon * (
                    np.power(self.sigma * this_distance, 12) - np.power(self.sigma * this_distance, 6))
        potential_energy = np.multiply(scaled, potential_energy)

        result = np.sum(potential_energy)

        del potential_energy
        del this_distance

        return result

    def calc_force(self, trajectory):
        """
        .. math:: -dU(r)/d(r) = [48\epsilon(\sigma^{12}/r^{13}) - 24\epsilon(\sigma^{6}/r^{7})]\cdot\overrightarrow{u}

        :param trajectory: coordinates

        :return: force matrix
        """
        # this_distance = self.pbcs_distance(trajectory)
        this_distance = distance.cdist(trajectory, trajectory)
        this_distance = np.where(this_distance == 0.0, 0.0, np.reciprocal(this_distance))

        vect = np.zeros((self.natoms, self.natoms, 3))
        for k in range(self.natoms):
            for j in range(self.natoms):
                tmp = trajectory[k] - trajectory[j]
                if k != j:
                    vect[k][j] = tmp / np.linalg.norm(tmp)
                else:
                    vect[k, j] = np.zeros(3)

        _force = (48 * self.epsilon * np.power(self.sigma, 12) * np.power(this_distance, 13)) - \
                 (24 * self.epsilon * np.power(self.sigma, 6) * np.power(this_distance, 7))

        force = _force.reshape(self.natoms, self.natoms, 1) * vect

        return force

    def calc_velocity(self, last_velocity, last_force, this_force):
        """
        Update velocity based on Verlet integration

        :param last_velocity: velocity of last step

        :param last_force: force on last step

        :param this_force: force on this step

        :return: velocity matrix
        """
        velocity = last_velocity + 0.5 * ((last_force + this_force) / self.mass) * self.deltaT
        return velocity

    def calc_position(self, last_position, last_velocity, last_force):
        """
        Update position based on Verlet integration

        :param last_position: position of last step

        :param last_velocity: velocity of last step

        :param last_force: force on this step

        :return: trajectory of this step
        """
        last_velocity = np.sum(last_velocity, axis=1)
        last_force = np.sum(last_force, axis=1)
        position = last_position + last_velocity * self.deltaT + 0.5 * (last_force / self.mass) * np.square(self.deltaT)
        position = np.where(position > self.system_size, position - self.system_size, position)
        return position

    def calc_kinetic_energy(self, this_velocity):
        """
        Calculate kinetic energies based on velocity and mass

        :param this_velocity: velocity of this step

        :return: kinetic energies matrix
        """
        kinetic = 0.5 * self.mass * np.square(np.sum(this_velocity, axis=1))
        kinetic = np.sum(kinetic)  # unit: kcal/mol
        return kinetic

    def run(self, steps):
        """
        Run simulations

        :param steps: simulations steps
        """
        step = 0
        while step < steps:
            if step == 0:
                last_position = self.topology
                last_velocity = np.random.uniform(-0.1, 0.1, size=(self.natoms, self.natoms, 3))
                self.kinetic_energies.append(self.calc_kinetic_energy(last_velocity))
                self.potential_energies.append(self.calc_potential_energy(last_position))
                this_force = self.calc_force(last_position)
                this_velocity = self.calc_velocity(last_velocity, this_force, this_force)
                last_force = this_force
                last_velocity = this_velocity
                self.velocities.append(last_velocity)
                self.trajectories.append(last_position)
                step += 1
            else:
                this_position = self.calc_position(last_position, last_velocity, last_force)
                self.potential_energies.append(self.calc_potential_energy(this_position))
                this_force = self.calc_force(this_position)
                this_velocity = self.calc_velocity(last_velocity, last_force, this_force)
                self.kinetic_energies.append(self.calc_kinetic_energy(this_velocity))
                last_force = this_force
                last_velocity = this_velocity
                last_position = this_position
                self.trajectories.append(last_position)
                self.velocities.append(last_velocity)
                step += 1
