import copy
import random
from statistics import mean, stdev

class GeneticAlgorithm(object):
    """
    A generic model of a genetic algorithm which can be overriden to specialise it
    towards a given task. See example.py for an example of this.

    ## Algorithm Parameters

    population_size: int, population size which is generated after mating
    mating_pool_size: int, population size of the fittest chromosomes selected for breeding

    # mating parameters
    mutation_coefficient: positive float, how much mutation should occur when transitioning generations
    mutation_step: absolute upper limit of how much each gene should be modified by
    crossover_strategy: what strategy to use when crossing over chromosomes.
                               success | generations until success (N=2000)
        INTERLEAVE_RATIO   = 0 : 0.989 | 11.012
        INTERLEAVE_RANDOM  = 1 : 0.971 | 11.479
        INTERLEAVE_FITNESS = 2 : 0.987 | 11.126
    crossover_survival: bool, keep highest fitness chromosomes across generations
    crossover_survival_limit: int, how many chromosomes to keep across generations

    # end conditions
    either one or both can be set

    generation_limit: int, how many generations until we exit
    convergence_similarity_limit: float, limit at which difference between max fitness'
                                         across subsequent generations is considered equal.
    convergence_similarity_occurrence_limit: int, after how many "equal" max fitness' across
                                                  generations we exit the spin.
    """

    INTERLEAVE_RATIO = 0
    INTERLEAVE_RANDOM = 1
    INTERLEAVE_FITNESS = 2

    def __init__(self, initial_population, debug=False):
        """
        :param initial_population: list of Chromosomes
        """
        self.population = initial_population
        self.generation = 0
        self.similarity_count = 0
        self.highest_fitness = 0
        self.fittest_candidate = None
        self.debug = debug

        # config parameters
        self.population_size = len(initial_population)
        self.mating_pool_size = self.population_size // 5

        self.mutation_coefficient = 1.0
        self.mutation_step = 1.0

        self.crossover_strategy = GeneticAlgorithm.INTERLEAVE_RATIO

        self.generation_limit = 100
        self.generation_fittest_survival_limit = 1

        self.convergence_similarity_limit = 0.1
        self.convergence_similarity_occurrence_limit = 2


    def fitness(self, chromosome):
        """
        Overwrite this with your own implementation
        :param chromosome: type Chromosome, to be assessed
        :return: fitness of chromosome
        """
        print("Please override the fitness function! Returning rand(0, 1) for now.")
        return random.randrange(0., 1.)

    def interleave(self, chromosome_a, chromosome_b, ratio=0.5, strategy=INTERLEAVE_RATIO):
        """
        :param chromosome_a: First Chromosome
        :param chromosome_b: Second Chromosome
        :param ratio: Point at which interleaving should happen, default 0.5 (1:1 ratio)
        :param strategy: either INTERLEAVE_RATIO, INTERLEAVE_RANDOM or INTERLEAVE_FITNESS
        :return: new Chromosome which is an offspring of chromosome_a and chromosome_b
        """
        if strategy == GeneticAlgorithm.INTERLEAVE_RATIO:
            return Chromosome(chromosome_a.split(start=0.0, end=ratio) + chromosome_b.split(start=ratio, end=1.0))
        elif strategy == GeneticAlgorithm.INTERLEAVE_RANDOM:
            return Chromosome([chromosome_a[i] if random.randint(0, 1) == 0 else chromosome_b[i]
                    for i in range(len(chromosome_a))])
        elif strategy == GeneticAlgorithm.INTERLEAVE_FITNESS:
            if not chromosome_a.fitness or not chromosome_b.fitness:
                print(chromosome_a, chromosome_b)
                raise ValueError("Chromosome does not have a fitness value!")
            return Chromosome([chromosome_a[i] if random.uniform(0, chromosome_a.fitness + chromosome_b.fitness)
                                       < chromosome_a.fitness else chromosome_b[i]
                    for i in range(len(chromosome_a))])

    def crossover(self, parent_pool):
        """
        Takes a pool of Chromosomes, and repeatedly interleaves 2 random Chromosomes' genes
        to produce children.
        :param parent_pool: list of Chromosome objects
        :return: list of Chromosome objects, of length population_size
        """
        children = []
        for index in range(self.population_size - self.generation_fittest_survival_limit):
            parent_a = parent_pool[random.randrange(0, len(parent_pool))]
            parent_b = parent_pool[random.randrange(0, len(parent_pool))]
            child = self.interleave(parent_a, parent_b, strategy=self.crossover_strategy)

            assert len(child) == len(parent_pool[0])
            children.append(child)

        return children

    def mutate(self, chromosome):
        """
        Mutates gene(s) in the input chromosome, mutation size and mutation count
        affected by mutation_step and mutation_coefficient
        :param chromosome: chromosome to be mutated, type Chromosome
        """
        chromosome[random.randrange(0, len(chromosome))] += random.uniform(-1., 1.) * self.mutation_step

    def mate(self, candidates):
        fittest_candidates = copy.deepcopy(candidates[:self.generation_fittest_survival_limit])
        if self.debug:
            print('Saving Fittest: ', fittest_candidates)
        mutations_remaining = self.mutation_coefficient * self.population_size
        if self.debug:
            print("CANDIDATES")
            print([(candidate, candidate.fitness) for candidate in candidates])
        while mutations_remaining > 0:
            self.mutate(candidates[random.randrange(0, len(candidates))])
            mutations_remaining -= 1
        if self.debug:
            print("AFTER MUTATION")
            print([(candidate, candidate.fitness) for candidate in candidates])
        next_population = self.crossover(candidates)
        if self.debug:
            print("AFTER CROSSOVER")
            print([(candidate, candidate.fitness) for candidate in next_population])
        return fittest_candidates + next_population

    def spin(self):
        """
        Spin until fittest candidate is found.
        Exits when either of the following is met:
            1) Generation limit is exceeded
            2) Similarity limit is exceeded
        :return: fittest candidate, type: Chromosome
        """

        if self.crossover_strategy not in [0, 1, 2]:
            self.crossover_strategy = GeneticAlgorithm.INTERLEAVE_RATIO

        while self.generation < self.generation_limit or not self.generation_limit:
            for chromosome in self.population:
                chromosome.fitness = self.fitness(chromosome)

            candidates = sorted(
                self.population,
                key=lambda candidate: candidate.fitness,
                reverse=True
            )[:self.mating_pool_size]

            fitnesses = [candidate.fitness for candidate in candidates]
            fitness_mean = mean(fitnesses)
            fitness_std = stdev(fitnesses)
            self.fittest_candidate = candidates[0]

            print("Generation %i: Highest fitness = %.2f, Mean fitness = %.2f (std=%.2f)"
                  % (self.generation, self.fittest_candidate.fitness, fitness_mean, fitness_std))
            print("Fittest Chromosome: " + str(self.fittest_candidate))

            if abs(self.fittest_candidate.fitness - self.highest_fitness) < self.convergence_similarity_limit:
                self.similarity_count += 1
            else:
                self.similarity_count = 0

            if self.highest_fitness < self.fittest_candidate.fitness:
                self.highest_fitness = self.fittest_candidate.fitness

            if self.similarity_count > self.convergence_similarity_occurrence_limit:
                print("Exiting due to similarity convergence", end="\n\n")
                return self.population

            #  set up next population
            self.population = self.mate(candidates)
            self.generation += 1
            print()

        return self.fittest_candidate


class Chromosome(list):
    def __init__(self, genes):
        super().__init__(genes)
        self.fitness = None

    def split(self, start=0.0, end=1.0):
        """
        :param start: float with value from 0 to 1.
        :param end: float with value from 0 to 1.
        :return: genome split at percentage determined in start and end
        """
        if start > end:
            raise ValueError("start > end: %.2f > %.2f" % (start, end))
        if start < 0.0 or start > 1.0 or end < 0.0 or end > 1.0:
            raise ValueError("start and/or end are out of bounds of 0.0 to 1.0: start=%.2f, end=%.2f" % (start, end))

        return self[int(len(self) * start):int(len(self) * end)]

    def __str__(self):
        return "[fitness={}, genome={}]".format(self.fitness, super().__str__())
