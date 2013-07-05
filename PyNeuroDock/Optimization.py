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
from time import time
import numpy as np
import pyopencl as cl
from pyopencl.clrandom import RanluxGenerator

DEBUG = False
DEBUG = True #bar

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
                #bar - start
                """
                individual.translation_gene = Axis3(2.056477, 5.846611, -7.245407)
                individual.rotation_gene = Quaternion(0.532211, 0.379383, 0.612442, 0.444674)
                torsions_gene_degrees = [-122.13, -179.41, \
                     -141.59,  177.29, \
                     -179.46,   -9.31, \
                     132.37,  -89.19, \
                     78.43,   22.22, \
                     71.37,   59.52]
                individual.torsions_gene = []
                for torsions_gene_degree in torsions_gene_degrees:
                    individual.torsions_gene.append((3.1415926535897931 / 180.0) * torsions_gene_degree)
                print individual.translation_gene
                print individual.rotation_gene
                print individual.torsions_gene
                """
                #bar - stop
                if dock.reset_pose(individual.translation_gene, \
                                   individual.rotation_gene, \
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

        # Create multiple population
        self.nomad = None
        self.settler = None

        # Define random number generator
        self.rng = LFSR(lfsr = 1070, bit_len = 48)

    def setup(self):
        self.lo_grid = self.dock.grid.field.lo
        self.hi_grid = self.dock.grid.field.hi
        self.ttl_torsions = self.dock.get_total_torsions()
        self.mutation_chance = 1.0 / self.ttl_torsions
        self.ttl_ligand_atoms = len(self.dock.ligand.ori_atoms)

        # Create multiple population
        self.nomad = self.Nomad()
        self.settler = self.Settler()

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
        pop_min_scores = []
        #return #bar  to check preprocessing time
        for community_idx in xrange(self.community_size):
            tic = time()
            # Nomad portion
            nomad_min_score = float("inf")
            self.nomad.setup(self.pop_size, self.ttl_torsions, \
                             self.lo_grid, self.hi_grid, \
                             self.rng)
            if DEBUG: print self.nomad
            for gen_idx in xrange(self.num_gen):
                mating_pool = self.select(self.nomad)
                self.nomad = self.reproduce(mating_pool, self.nomad)
            nomad_min_score = self.nomad.scores.minimum()

            # Settler portion
            settler_min_score = float("inf")
            self.settler.individuals = deepcopy(self.nomad.individuals)
            if DEBUG: print self.settler
            for gen_idx in xrange(self.num_gen):
                mating_pool = self.select(self.settler)
                self.settler = self.reproduce(mating_pool, self.settler)
            if DEBUG: print self.settler
            settler_min_score = self.settler.scores.minimum()

            pop_min_scores.append([nomad_min_score, settler_min_score])
            if DEBUG:
                print "Current Minimum Scores: %f, %f" % (nomad_min_score, \
                                                          settler_min_score)
                print "Community Minimum Scores: %s" % pop_min_scores
            toc = time()
            print "Elapsed time community %3d: %10.2f" % (community_idx + 1, \
                                                          toc - tic)

class GeneticAlgorithmOpenCL(GeneticAlgorithm):
    class Population:
        def __init__(self, size = 0):
            self.size = size
            # Matrix of i by j for individuals and genes (DNA)
            self.individuals = None
            self.individuals_buf = None
            self.scores = GeneticAlgorithmOpenCL.Scores()

        def __repr__(self):
            self.individuals = self.individuals_buf.get()

            ret = "Individuals:\n"
            for idx, individual in enumerate(self.individuals):
                ret += "[%3d] " % (idx + 1)
                ret += "%6.2f %6.2f %6.2f | " % (individual[0], \
                                                 individual[1], \
                                                 individual[2])
                ret += "%6.2f %6.2fi %6.2fj %6.2fk | " % (individual[3], \
                                                          individual[4], \
                                                          individual[5], \
                                                          individual[6])
                for torsion in individual[7:]:
                    ret += " %5.2f" % torsion
                ret += "\n"
            return ret

        def setup(self, ttl_torsions = 0, dock = None, \
                  cl_queue = None, rng = None, cl_prg = None):
            dna_size = 3 + 4 + int(ttl_torsions)
            self.individuals_buf = cl.array.zeros(cl_queue, \
                                                  (self.size, dna_size), \
                                                  dtype = float)
            rng.fill_uniform(self.individuals_buf)
            # Construct individuals
            cl_prg.construct_individuals(cl_queue, (self.size,), None, \
                                         self.individuals_buf.data, \
                                         dock.lo_grid_buf, dock.dist_grid_buf, \
                                         dock.ttl_torsions_buf)

        def scoring(self, dock = None, \
                    cl_ctx = None, cl_queue = None):
            dock.reset_poses(self.size, self.individuals_buf, \
                             cl_ctx, cl_queue)
            dock.calc_energy()

        def crossover(self, parents_idx, ttl_torsions, rng):
            return None

        def mutate(self, individual, mutation_chance, \
                   lo_grid, hi_grid, ttl_torsions, rng):
            return None

    class Settler(Population):
        def __init__(self, size = 0):
            GeneticAlgorithmOpenCL.Population.__init__(self, size)

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
        def __init__(self, size = 0):
            GeneticAlgorithmOpenCL.Population.__init__(self, size)

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

    def __init__(self, dock = None):
        GeneticAlgorithm.__init__(self, dock)
        self.dock = dock
        self.community_size = 0         # Community size
        self.pop_size = 0               # Population size
        self.num_gen = 0                # Number of generations
        self.max_inherited_prob = 12    # Maximum inhereted probability
        self.mutation_chance = 0.0
        self.ttl_ligand_atoms = 0

        # Create multiple population
        self.nomad = None
        self.settler = None

        # OpenCL
        self.cl_ctx = None
        self.cl_queue = None
        self.rng = None
        self.cl_filename = "./OpenCL/GeneticAlgorithm.cl"
        self.cl_prg = None
        # OpenCL buffer
        self.max_inherited_prob_np = np.array([], dtype = int)
        self.max_inherited_prob_buf = None
        self.normalizer_np = np.array([], dtype = int)
        self.normalizer_buf = None
        self.chances_np = None
        self.chances_buf = None

    def setup(self):
        # Create multiple population
        self.nomad = self.Nomad(self.pop_size)
        self.settler = self.Settler(self.pop_size)
        
        self.ttl_torsions = self.dock.get_total_torsions()
        self.mutation_chance = 1.0 / self.ttl_torsions
        self.ttl_ligand_atoms = len(self.dock.ligand.ori_atoms)

        # OpenCL setup
        self.cl_ctx = cl.Context(dev_type = cl.device_type.GPU)
        self.cl_queue = cl.CommandQueue(self.cl_ctx)
        self.rng = RanluxGenerator(self.cl_queue)
        # Read OpenCL code
        fh = open(self.cl_filename, 'r')
        cl_code = "".join(fh.readlines())
        self.cl_prg = cl.Program(self.cl_ctx, cl_code).build()
        self.dock.setup_opencl(self.cl_ctx, self.cl_queue)

        # Setup OpenCL device buffer
        mf = cl.mem_flags
        self.max_inherited_prob_np = np.array([self.max_inherited_prob], \
                                              dtype = int)
        self.max_inherited_prob_buf = cl.Buffer(self.cl_ctx, \
                                                mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                                hostbuf = self.max_inherited_prob_np)
        self.normalizer_np = np.array([self.ttl_ligand_atoms], dtype = int)
        self.normalizer_buf = cl.Buffer(self.cl_ctx, \
                                        mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                        hostbuf = self.normalizer_np)
        self.chances_buf = cl.array.zeros(self.cl_queue, (self.pop_size), \
                                          dtype = int)
        self.dock.setup_opencl_buffer(self.pop_size, \
                                      self.cl_ctx, self.cl_queue)

    def select(self, population):
        # Get individual scores
        population.scoring(self.dock, self.cl_ctx, self.cl_queue)
        self.cl_prg.calc_chances(self.cl_queue, (population.size,), None, \
                                 self.dock.e_totals_buf.data,
                                 self.normalizer_buf,
                                 self.max_inherited_prob_buf,
                                 self.chances_buf.data)
        if DEBUG:
            self.chances_np = self.chances_buf.get()
            print self.chances_np

        mating_pool = []
        for idx, chance in enumerate(self.chances_np):
            for i in xrange(chance):
                mating_pool.append(idx)

        if DEBUG:
            print mating_pool
            print len(mating_pool)

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
        pop_min_scores = []
        #return #bar to check preprocessing itme
        for community_idx in xrange(self.community_size):
            tic = time()
            # Nomad portion
            nomad_min_score = float("inf")
            self.nomad.setup(self.ttl_torsions, self.dock, \
                             self.cl_queue, self.rng, self.cl_prg)
            if DEBUG: print self.nomad
            for gen_idx in xrange(self.num_gen):
                mating_pool = self.select(self.nomad)
                self.nomad = self.reproduce(mating_pool, self.nomad)
            nomad_min_score = self.nomad.scores.minimum()

            # Settler portion
            settler_min_score = float("inf")
            settler.individuals = deepcopy(self.nomad.individuals)
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
            toc = time()
            print "Elapsed time community %3d: %10.2f" % (community_idx + 1, \
                                                          toc - tic)
