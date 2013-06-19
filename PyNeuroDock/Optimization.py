# Copyright (C) 2013 by Eka A. Kurniawan
# eka.a.kurniawan(ta)gmail(tod)com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# References:
#  - AutoDock 4.2.3 Source Code (call_glss.cc)
#    http://autodock.scripps.edu
#  - The Nature of Code by Daniel Shiffman (Chapter 9. The Evolution of Code)
#    http://natureofcode.com/book/chapter-9-the-evolution-of-code/

from Axis3 import Axis3
from Quaternion import Quaternion
from Constants import DEG2RAD
from LFSR import LFSR
from math import log
from copy import deepcopy

DEBUG = True

class GeneticAlgorithm:
    class Individual:
        def __init__(self, ttl_torsions = 0, \
                     lo_grid = Axis3(), hi_grid = Axis3(), \
                     rng = None):
            self.translation_gene = Axis3()
            self.rotation_gene = Quaternion()
            self.torsions_gene = []
            
            self.random_translation(lo_grid, hi_grid, rng)
            self.random_rotation(rng)
            self.random_torsions(ttl_torsions, rng)

        def random_translation(self, lo_grid, hi_grid, rng):
            lo_x, lo_y, lo_z = lo_grid.xyz
            hi_x, hi_y, hi_z = hi_grid.xyz
            self.translation_gene.x = lo_x + (rng.zero_to_one() * (hi_x - lo_x))
            self.translation_gene.y = lo_y + (rng.zero_to_one() * (hi_y - lo_y))
            self.translation_gene.z = lo_z + (rng.zero_to_one() * (hi_z - lo_z))

        def random_rotation(self, rng):
            self.rotation_gene.uniform(rng)

        def random_torsions(self, ttl_torsions, rng):
            self.torsions_gene = []
            for i in xrange(ttl_torsions):
                self.torsions_gene.append(rng.neg_pi_to_pi())

    class Population:
        def __init__(self):
            self.individuals = []
            self.scores = GeneticAlgorithm.Scores()

        def __repr__(self):
            ret = "Individuals:\n"
            for individual in self.individuals:
                ret += "%s  %s |" % (individual.translation_gene, \
                                     individual.rotation_gene)
                for torsion in individual.torsions_gene:
                    ret += " %5.2f" % torsion
                ret += "\n"
            return ret

        def setup(self, size = 0, ttl_torsions = 0, \
                  lo_grid = Axis3(), hi_grid = Axis3(), \
                  rng = None):
            self.individuals = []
            for i in xrange(size):
                individual = GeneticAlgorithm.Individual(ttl_torsions, \
                                                         lo_grid, hi_grid, \
                                                         rng)
                self.individuals.append(individual)

        def scoring(self, dock = None):
            self.scores = GeneticAlgorithm.Scores()
            for individual in self.individuals:
                if dock.reset_pose(individual.translation_gene, \
                                   individual.rotation_gene,
                                   individual.torsions_gene):
                    self.scores.append(dock.calc_energy())
                else:
                    self.scores.append(float("inf"))
            return self.scores

        def crossover(self, parents_idx, ttl_torsions, rng):
            return None

        def mutate(self, individual, mutation_chance, \
                   lo_grid, hi_grid, ttl_torsions, rng):
            return None

    class Settler(Population):
        def __init__(self):
            GeneticAlgorithm.Population.__init__(self)

        def crossover(self, parents_idx, ttl_torsions, rng):
            p1_idx = parents_idx[0]
            p2_idx = parents_idx[1]
            individual = deepcopy(self.individuals[p1_idx])
            # Crossing over translation gene
            if rng.zero_to_one() > 0.5:
                individual.translation_gene = deepcopy(self.individuals[p2_idx].translation_gene)
            # Crossing over rotation gene
            if rng.zero_to_one() > 0.5:
                individual.rotation_gene = deepcopy(self.individuals[p2_idx].rotation_gene)
            # Crossing over torsions gene
            for i in xrange(ttl_torsions):
                if rng.zero_to_one() > 0.5:
                    individual.torsions_gene[i] = self.individuals[p2_idx].torsions_gene[i]
            return individual

        def mutate(self, individual, mutation_chance, \
                   lo_grid, hi_grid, ttl_torsions, rng):
            if rng.zero_to_one() < mutation_chance:
                # Mutating translation gene
                if rng.zero_to_one() < 0.25:
                    individual.random_translation(lo_grid, hi_grid, rng)
                # Mutating rotation gene
                if rng.zero_to_one() < 0.25:
                    individual.random_rotation(rng)
                # Mutating torsions gene
                for i in xrange(ttl_torsions):
                    if rng.zero_to_one() < 0.25:
                        individual.torsions_gene[i] = rng.neg_pi_to_pi()
            return individual

    class Nomad(Population):
        def __init__(self):
            GeneticAlgorithm.Population.__init__(self)

        def crossover(self, parents_idx, ttl_torsions, rng):
            p1_idx = parents_idx[0]
            p2_idx = parents_idx[1]
            individual = deepcopy(self.individuals[p1_idx])
            # Crossing over translation gene
            if rng.zero_to_one() > 0.5:
                individual.translation_gene.x = self.individuals[p2_idx].translation_gene.x
            if rng.zero_to_one() > 0.5:
                individual.translation_gene.y = self.individuals[p2_idx].translation_gene.y
            if rng.zero_to_one() > 0.5:
                individual.translation_gene.z = self.individuals[p2_idx].translation_gene.z
            # Crossing over rotation gene
            if rng.zero_to_one() > 0.5:
                individual.rotation_gene.a = self.individuals[p2_idx].rotation_gene.a
            if rng.zero_to_one() > 0.5:
                individual.rotation_gene.b = self.individuals[p2_idx].rotation_gene.b
            if rng.zero_to_one() > 0.5:
                individual.rotation_gene.c = self.individuals[p2_idx].rotation_gene.c
            if rng.zero_to_one() > 0.5:
                individual.rotation_gene.d = self.individuals[p2_idx].rotation_gene.d
            # Crossing over torsions gene
            for i in xrange(ttl_torsions):
                if rng.zero_to_one() > 0.5:
                    individual.torsions_gene[i] = self.individuals[p2_idx].torsions_gene[i]
            return individual

        def mutate(self, individual, mutation_chance, \
                   lo_grid, hi_grid, ttl_torsions, rng):
            if rng.zero_to_one() < mutation_chance:
                # Mutating translation gene
                if rng.zero_to_one() < 0.75:
                    individual.random_translation(lo_grid, hi_grid, rng)
                # Mutating rotation gene
                if rng.zero_to_one() < 0.75:
                    individual.random_rotation(rng)
                # Mutating torsions gene
                for i in xrange(ttl_torsions):
                    if rng.zero_to_one() < 0.75:
                        individual.torsions_gene[i] = rng.neg_pi_to_pi()
            return individual

    class Scores(list):
        def __init__(self, *args):
            list.__init__(self, *args)

        def minimum(self):
            return min(self)

        # Normalized scores to the total ligand atoms
        def normalize(self, normalizer):
            normalized_scores = []
            for score in self:
                normalized_scores.append(float(score) / normalizer)
            return normalized_scores

    def __init__(self, dock = None):
        self.dock = dock
        self.lo_grid = Axis3()
        self.hi_grid = Axis3()
        self.ttl_torsions = 0
        self.ttl_ligand_atoms = 0

        self.community_size = 0         # Community size
        self.pop_size = 0               # Population size
        self.num_gen = 0                # Number of generations
        self.max_inherited_prob = 12    # Maximum inhereted probability

        # Define random number generator
        self.rng = LFSR(lfsr = 1070, bit_len = 48)

    def setup(self):
        self.lo_grid = self.dock.grid.field.lo
        self.hi_grid = self.dock.grid.field.hi
        self.ttl_torsions = self.dock.get_total_torsions()
        self.mutation_chance = 1.0 / self.ttl_torsions
        self.ttl_ligand_atoms = len(self.dock.ligand.ori_atoms)

    # Create multiple population with predefined number of individuals
    def initialize(self):
        nomad = self.Nomad()
        settler = self.Settler()
        return nomad, settler

    def select(self, population):
        # Get individual scores
        scores = population.scoring(self.dock)
        # Create mating pool from the scores
        mating_pool = []
        for idx, score in enumerate(scores.normalize(self.ttl_ligand_atoms)):
            # Use probabilistic method to select individual into mating pool
            if score < 0:
                chances = self.max_inherited_prob
            else:
                power = log(abs(score))
                if power < self.max_inherited_prob:
                    chances = int(self.max_inherited_prob - power)
                else:
                    chances = 1
            # Fill in the mating pool
            for i in xrange(chances):
                mating_pool.append(idx)
        return mating_pool

    def pick_parents(self, mating_pool):
        parents = []
        for i in xrange(2):
            parents.append(mating_pool[int(self.rng.zero_to_one() * len(mating_pool))])
        return parents

    def reproduce(self, mating_pool, population):
        new_population = deepcopy(population)
        new_population.individuals = []
        for i in xrange(self.pop_size):
            parents_idx = self.pick_parents(mating_pool)
            individual = population.crossover(parents_idx, self.ttl_torsions, self.rng)
            individual = population.mutate(individual, self.mutation_chance, \
                                           self.lo_grid, self.hi_grid, \
                                           self.ttl_torsions, self.rng)
            new_population.individuals.append(individual)
        return new_population

    def run(self):
        self.setup()
        self.initialize()
        nomad, settler = self.initialize()
        pop_min_scores = []
        for community_idx in xrange(self.community_size):
            # Nomad portion
            nomad_min_score = float("inf")
            nomad.setup(self.pop_size, self.ttl_torsions, \
                        self.lo_grid, self.hi_grid, \
                        self.rng)
            if DEBUG: print nomad
            for gen_idx in xrange(self.num_gen):
                mating_pool = self.select(nomad)
                nomad = self.reproduce(mating_pool, nomad)
            nomad_min_score = nomad.scores.minimum()

            # Settler portion
            settler_min_score = float("inf")
            settler.individuals = deepcopy(nomad.individuals)
            if DEBUG: print settler
            for gen_idx in xrange(self.num_gen):
                mating_pool = self.select(settler)
                settler = self.reproduce(mating_pool, settler)
            if DEBUG: print settler
            settler_min_score = settler.scores.minimum()

            pop_min_scores.append([nomad_min_score, settler_min_score])
            if DEBUG:
                print "Current Minimum Scores: %f, %f" % (nomad_min_score, \
                                                          settler_min_score)
                print "Community Minimum Scores: %s" % pop_min_scores

