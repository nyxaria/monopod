#!/usr/bin/env python3

import random
import string

from model import GeneticAlgorithm, Chromosome
from shapely.geometry import Polygon
import shapely.affinity as shapely_affinity


class FourBarLinkageSolver(GeneticAlgorithm):
    """
    An implemenation of a genetic algorithm specialised to generate a 4 bar linkage
    for a given input end effector path.

    ___ Documentation ___

    In a chromosome, a 4 bar linkage is defined as follows: [l1, l2, l3].
    Each variable is a positive, non-zero float

                                   ~
                                  / l5
                             l2  /
                           _____/
    example:          l1  /      \  l3
                         /        \
                        *----------
                             l4

    * represents the rotational joint, about which l1 rotates.
    ~ represents the end effector, which follows the target geometry.
    l4 is the fixed linkage. We do not define it in the chromosome as it is implied by l1-3.

    The target geometry is defined the variable target_geometry, passed in the init function.
    This has the following format: [(x0,y0), (x1,y1), ... , (xn, yn)]
    The first and last points will be stitched to form a closed polygon.
    Furthermore, the variable target_geometry_offset is the translation between the rotational
    joint to the centroid of the polygon formed by the target_geometry polygon.

    """

    def __init__(self, initial_population, target_geometry, target_geometry_offset, debug=False):
        super().__init__(initial_population, debug=debug)
        if(type(target_geometry) is not list or
           type(target_geometry[0]) is not tuple or
           len(target_geometry[0]) != 2 or
           type(target_geometry[0][0]) is not float):
            raise ValueError("target_geometry must be of the following format: [(x0,y0), (x1,y1), ... , (xn, yn)]" +
                             "All values are floats.")

        if((type(target_geometry_offset) is not list and type(target_geometry) is not tuple) or
           len(target_geometry_offset) != 2 or
           type(target_geometry[0]) is not float):
            raise ValueError("target_geometry_offset must be of the following format: [x_off, y_off]" +
                             "All values are floats. Can also be a tuple.")
        self.target_geometry = shapely_affinity.translate(Polygon(target_geometry),
                                                          xoff=target_geometry_offset[0],
                                                          yoff=target_geometry_offset[1])

    def fitness(self, chromosome):  # Fitness = How closely the chromosome matches the target word
        # inital checks to verify that this proposed 4 bar linkage is valid
        # l1 + l4 < l2 + l3
        if(chromosome[0] + chromosome)

        # TODO generate polygon from the chromosome

        chromosome_geom = self.get_geometry_from_chromosome(chromosome)

        Polygon([(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)])

        A.contains(B)

        A.difference(B).wkt
        'POLYGON ((60 350, 430 350, 430 30, 60 30, 60 350), (229 255, 156 184, 130 60, 230 110, 229 255))'
        B.difference(A).wkt
        'GEOMETRYCOLLECTION EMPTY'
        mapping(A.difference(B))
        {'type': 'Polygon', 'coordinates': (((60.0, 350.0), (430.0, 350.0), (430.0, 30.0), (60.0, 30.0), (60.0, 350.0)), ((
            229.0, 255.0), (156.0, 184.0), (130.0, 60.0), (230.0, 110.0), (229.0, 255.0)))}
        mapping(B.difference(A))
        {'type': 'GeometryCollection', 'geometries': []}

        if(area_overlap == 0):  # win condition; matches perfectly
            return 10**10000
        return 1/area_overlap  # highest value = fitness, so we inverse the difference

    def mutate(self, chromosome):  # Override default implementation to never get negative/zero lengths
        index = random.randrange(0, len(chromosome))
        chromosome[index] = abs(
            chromosome[index] + random.randint(-1, 1) * self.mutation_step)

        # check for zero length
        if chromosome[index] < 0.0001:
            chromosome[index] = random.randint(0, 1) * self.mutation_step


def main():
    population_size = 10000  # the higher the pop size, the less generations
    target_word = "potatoes"

    print("\nTarget Word: %s\n" % target_word)
    initial_population = [Chromosome([random.choice(string.ascii_lowercase) for i in range(len(target_word))])
                          for y in range(population_size)]  # create Chromosomes and their Genes, randomly

    ga = WordSolverGeneticAlgorithm(
        initial_population, target_word, debug=False)

    # The coefficients can be altered here

    # ga.mutation_coefficient = 0.5
    # ga.mutation_step = 1
    # ga.mating_pool_size = population_size//10

    ga.spin()

    print("Done. Fittest Candidate:")
    print(ga.fittest_candidate, "fitness = " +
          str(ga.fittest_candidate.fitness))


if __name__ == "__main__":
    main()
