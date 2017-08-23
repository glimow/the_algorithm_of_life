import random
import numpy as np
random.seed()


#This parameter controls average lifetime of individuals
LIFETIME_FACTOR = 10
#This parameter controls max and default hunger
HUNGER_FACTOR = 10



def mutate_gene(gene, mutation_threshold):
    seed = random.random()
    if seed < mutation_threshold:
        seed = seed*2 - mutation_threshold # seed must be btw -1 and 1
        # print seed
        return max(min(gene + seed,1),-1)
    else:
        return gene
# mutate_vector = np.vectorize(mutate_gene)

class Individual:

    # Each allele is represented by a number between -1 and 1
    #genotype = []


    genes_definition = [
    ["carnivorous",	"herbivorous"],
    ["speed",	"resistance"],
    ["flee",	"battle"],
    ["eat",	"reproduce"],
    ["n° of children",	"lifetime"],
    ["combat",	"hunger"],]

    genotype_size = len(genes_definition)

    mutation_treshold = 0.9
    mutation_treshold = 0.2

    def __init__(self, genotype):
        self.genotype = np.array(genotype)
        self.life = self.positive_gene(1)
        self.lifetime = int(self.positive_gene(4)*LIFETIME_FACTOR)
        self.age = 0
        self.hunger = self.positive_gene(5)*HUNGER_FACTOR


    def distance(self, other):
        return np.linalg.norm(self.genotype - other.genotype)

    def eat_other(self, other):
        self.hunger = min( self.hunger + other.lifetime*self.negative_gene(0), self.positive_gene(5)*HUNGER_FACTOR)

    def positive_gene(self, gene):
        return (self.genotype[gene]+1)/2*10

    def negative_gene(self, gene):
        return (-self.genotype[gene]+1)/2*10

    # return an array version of the Individual
    def to_array(self):
        return self.genotype

    # mutate the current individual
    def mutate(self):
        new_genotype = np.array([])
        # For each gene, we try to mutate this gene. Probably not the best way
        # to do it though
        for gene in self.genotype:
            mutation = mutate_gene(gene, self.mutation_treshold)
            # print mutation
            new_genotype = np.append(new_genotype,mutation)
        # We update the genotype of this individual
        self.genotype = new_genotype

    def to_genotype(self):
        return self.genotype

    # make children !
    def crossover(self, bae):
        parent1 = self.to_genotype()
        parent2 = bae.to_genotype()
        child = []
        # for each gene, launch a roll to chose from which parent to take it
        for gene1, gene2 in zip(parent1,parent2):
            seed = random.random()
            if seed > 0.5:
                child.append(gene1)
            else:
                child.append(gene2)
        child = Individual(child)
        #some random mutation
        child.mutate()
        return child

    def __str__(self):
        return self.genotype.__str__()

    def __repr__(self):
        return self.genotype.__str__()

if __name__ == "__main__":

    patient_0 = Individual(np.random.random_sample(Individual.genotype_size)*2-1)
    patient_1 = Individual(np.random.random_sample(Individual.genotype_size)*2-1)
    print("################# TESTING GENETICS ################")
    print("parent1", patient_0)
    patient_0.mutate()
    print("parent1, mutated", patient_0)
    print("parent2", patient_1)
    children = patient_0.crossover(patient_1)
    print("children", children)
    print("distance", patient_0.distance(patient_1))
    print("distance_children", patient_0.distance(children))
    print("################# TESTING ACTIONS ################")
    patient_0.hunger = patient_0.hunger/2.0
    print("0's hunger",patient_0.hunger )
    patient_0.eat_other(patient_1)
    print("0 eat 1", patient_0.hunger)
