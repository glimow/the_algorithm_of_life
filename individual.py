import random
import numpy as np
random.seed()

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
    genotype_size = 4
    mutation_treshold = 0.9
    mutation_treshold = 0.1

    def __init__(self, genotype):
        self.genotype = np.array(genotype)

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
            print mutation
            new_genotype = np.append(new_genotype,mutation)
        # We update the genotype of this individual
        self.genotype = new_genotype

    def __str__(self):
        return self.genotype.__str__()

if __name__ == "__main__":

    patient_0 = Individual([0.4,0.2,-0.4,0.3])
    print(patient_0)
    patient_0.mutate()
    print(patient_0)
