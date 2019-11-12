#!/usr/bin/env python3

import random
import string

from model import GeneticAlgorithm, Chromosome


class WordSolverGeneticAlgorithm(GeneticAlgorithm):
    """
    An implemenation of a genetic algorithm specialised to converge towards a
    target word.
    """
    def __init__(self, initial_population, target_word, debug=False):
        super().__init__(initial_population, debug=debug)
        self.target_word = target_word

    def fitness(self, chromosome):  # Fitness = How closely the chromosome matches the target word
        count = 0

        for index in range(len(self.target_word)):
            if self.target_word[index] == chromosome[index]:
                count += 1
        return count

    def mutate(self, chromosome):  # Override default implementation as we handle chars, not floats, as genes
        index = random.randrange(0, len(chromosome))
        chromosome[index] = chr(ord(chromosome[index]) + random.randint(-1, 1) * int(self.mutation_step))


def main():
    population_size = 10000  # the higher the pop size, the less generations
    target_word = "potatoes"

    print("\nTarget Word: %s\n" % target_word)
    initial_population = [Chromosome([random.choice(string.ascii_lowercase) for i in range(len(target_word))])
                          for y in range(population_size)]  # create Chromosomes and their Genes, randomly

    ga = WordSolverGeneticAlgorithm(initial_population, target_word, debug=False)

    # The coefficients can be altered here

    # ga.mutation_coefficient = 0.5
    # ga.mutation_step = 1
    # ga.mating_pool_size = population_size//10

    ga.spin()

    print("Done. Fittest Candidate:")
    print(ga.fittest_candidate, "fitness = " + str(ga.fittest_candidate.fitness))


if __name__ == "__main__":
    main()
