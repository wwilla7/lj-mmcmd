#!/usr/bin/env python
"""
This module is used for carrying out a simple Metropolis Monte Carlo simulation of Lennard Jones particles
"""
import numpy as np
import pint
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity


class MCLJ:
    def __init__(self):
        self.R = Q_(0.001985875, "kcal/mol/K").magnitude
        self.p_accept = 1
        self.accept = False
        self.trajectories = []
        self.potential_energies = []

    @property
    def epsilon(self):
        """
        Argon: Q_(0.238, "kcal/mol").magnitude
        :return:
        """
        try:
            return self._epsilon
        except AttributeError:
            print("Set epsilon!")

    @epsilon.setter
    def epsilon(self, value):
        if isinstance(value, pint.Quantity):
            self._epsilon = value.magnitude
        elif isinstance(value, float):
            self._epsilon = Q_(value, "kcal/mol").magnitude
            print("Default unit of epsilon, kcal/mol")
        else:
            print("Invalid epsilon", type(value))

    @property
    def sigma(self):
        """
        Argon: Q_(3.405, ureg.angstrom).magnitude
        :return:
        """
        try:
            return self._sigma
        except AttributeError:
            print("Set sigma!")

    @sigma.setter
    def sigma(self, value):
        if isinstance(value, pint.Quantity):
            self._sigma = value.magnitude
        elif isinstance(value, float):
            self._sigma = Q_(value, ureg.angstrom).magnitude
            print("Default unit of sigma, angstrom")
        else:
            print("Invalid sigma", type(value))


    @property
    def temperature(self):
        """
        Default 298 Kelvin
        :return:
        """
        try:
            return self._temperature
        except AttributeError:
            print("Set Temperature!")

    @temperature.setter
    def temperature(self, value):
        if isinstance(value, pint.Quantity):
            self._temperature = value.magnitude
        elif isinstance(value, (int, float)):
            self._temperature = Q_(value, ureg.kelvin).magnitude
        else:
            print("Invalid temperature", type(value))

    @property
    def system_size(self):
        """
        Unit: angstrom
        :return:
        """
        try:
            return self._system_size
        except AttributeError:
            print("Set system sizei!")

    @system_size.setter
    def system_size(self, value):
        self._system_size = value

    @property
    def nparticles(self):
        """
        Number of particles in this system
        :return:
        """
        try:
            return self._nparticles
        except AttributeError:
            print("Define number of particles in the system")

    @nparticles.setter
    def nparticles(self, value):
        self._nparticles = value

    def new_positions(self):
        """
        Generate a new position
        :return:
        """
        return np.random.uniform(0, self.system_size, size=(self.nparticles, 3))
        # return np.random.uniform(-self.system_size, self.system_size, size=(self.nparticles, 3))

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
        :param trajectory: coordinates
        :return: distance matrix with periodic boundary conditions
        """
        this_distance = np.zeros((self.nparticles, self.nparticles))
        for atom1 in range(self.nparticles):
            for atom2 in range(self.nparticles):
                this_distance[atom1, atom2] = self._pbcs(trajectory[atom1], trajectory[atom2])
        this_distance = np.tril(this_distance)
        return this_distance

    def calc_potential_energy(self, trajectory):
        """
        calculate the potential energy of lennard jones fluids
        V(r) = 4\epsilon[(\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6]
        switching function: S = 1 - 6x^5 + 15x^4 - 10x^3
        http://docs.openmm.org/latest/userguide/theory.html#lennard-jones-interaction
        :param trajectory: coordinates
        :return:
        """
        this_distance = self.pbcs_distance(trajectory)
        shift = lambda x: 1 - 6 * np.power(x, 5) + 15 * np.power(x, 4) - 10 * np.power(x, 3)
        this_distance = np.where(this_distance > 3*self.sigma, 0.0, this_distance)
        scaled = np.where(this_distance > 2*self.sigma, shift(((this_distance - 2*self.sigma) / (self.sigma))), 1.0)
        this_distance = np.where(this_distance == 0.0, 0.0, np.reciprocal(this_distance))

        potential_energy = 4 * self.epsilon * (np.power(self.sigma * this_distance, 12) - np.power(self.sigma * this_distance, 6))
        potential_energy = np.multiply(scaled, potential_energy)

        result = np.sum(potential_energy)

        del potential_energy
        del this_distance

        return result

    def _possibility(self, pot_energy1, pot_energy2):
        """
        P_{accept} = \frac{P_{trial}}{P_i}
                   = \frac{e^{-U_{trial}/RT}}{e^{-U/RT}}
                   = e^{-(U_{trial} - U_i)/RT}
        :param pot_energy1: Potential Energies of last position
        :param pot_energy2: Potential Energies of this position
        :return:
        """
        #possibility = np.power(np.e, -(pot_energy2 - pot_energy1) / (self.R * self.temperature))
        possibility = np.exp(-(pot_energy2 - pot_energy1)/(self.R*self.temperature))

        return possibility

    def decision_maker(self, pot_energy1, pot_energy2):
        """
        use Metropolis Monte Carlo algorithm to decide if the trial run will be accepted
        :param pot_energy1: Potential Energies of last position
        :param pot_energy2: Potential Energies of this position
        :return:
        """
        if pot_energy1 > pot_energy2:
            self.accept = True
        else:
            p_accept = self._possibility(pot_energy1, pot_energy2)
            mc_nu = np.random.uniform(0, 1, size=(1))[0]
            if mc_nu < p_accept:
                self.accept = True
            else:
                self.accept = False
