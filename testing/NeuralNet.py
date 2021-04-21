import numpy as np
from settings import *
from copy import copy
import json, codecs, random, logging
from helpers import *
from os import path

logging.basicConfig(filename="NeuralNetwork.log", level=logging.DEBUG)

class Neural_network():
    def __init__(self):
        self.input_size = 24
        self.hidden_layer_size = 16
        self.output_layer_size = 4
        self.genetics_dir = path.join(GAME_FOLDER, "Genetics")
        self.score_dir = path.join(GAME_FOLDER, "Scores")
        self.reference_file = path.join(GAME_FOLDER, "reference.json")

        self.mutation_rate = MUTATION_RATE
        self.crossover_rate = CROSSOVER_RATE
        self.genome_size = GENOME_SIZE
        self.selection_num = int(self.genome_size/3)
        self.activation_func = ACTIVATION_FUNCTION
        self.cloning_enabled = CLONING_ENABLED
        self.learning_rate = LEARNING_RATE

    def normalize(self, raw_inputs):
        normalized_input = []
        for x in raw_inputs:
            normalized_input.append(x/TILESIZE)

        return normalized_input

    def forwardProp(self, X_inputs, chromosome):
        '''
        Runs through forward propagation math and returns a list of the outputs of the
        neural network.

        Parameters: (list X_inputs, dict chromosome)
        '''

        W1 = chromosome["W1"]
        W2 = chromosome["W2"]
        X = None

        if NORMALIZE_INPUTS:
            X = self.normalize(X_inputs)
        elif not NORMALIZE_INPUTS:
            X = X_inputs

        Z2 = np.dot(X, W1) * self.learning_rate                     #Formula 1
        a2 = ''
        if self.activation_func == "SIGMOID":
            a2 = sigmoid(Z2)                       #Formula 2
        elif self.activation_func == "TANH":
            a2 = tanhActivation(Z2)                #Formula 2
        Z3 = np.dot(a2, W2) * self.learning_rate                    #Formula 3
        y_hat = ''
        if self.activation_func == "SIGMOID":
            y_hat = sigmoid(Z3)                    #Formula 4
        elif self.activation_func == "TANH":
            y_hat = tanhActivation(Z3)                    #Formula 4

        return y_hat

    def selection_elitism(self):
        current_species = self.get_current_species()
        current_genome = self.get_current_genome()
        genome = self.readFile(path.join(self.score_dir, "s_"+str(current_species)+"_genome_"+str(current_genome)+".json"))
        score_data = genome["scores"].copy()

        final_selected = []
        for i in range(self.selection_num):
            max_num = max(score_data)
            max_num_idx = score_data.index(max_num)
            final_selected.append([max_num_idx, max_num])
            score_data.pop(max_num_idx)
        #return a list of the selected chromosomes
        return final_selected

    def selection_Rchance(self):
        '''
        Select the fittest chromosomes in a genome and return a list of
        chromosome pairs. Finds and returns the top 4 highest chromosomes in
        each genome in order from best to worst.
        Parameters: (dict genome)
        '''
        current_species = self.get_current_species()
        current_genome = self.get_current_genome()
        genome = self.readFile(path.join(self.score_dir, "s_"+str(current_species)+"_genome_"+str(current_genome)+".json"))
        score_data = genome["scores"]


        #convert score_data list to dictionary with corresponding chromosome number as keys
        #Also find the total sum of all of the scores in the genome
        score_sum = 0
        chromosome_dict = []
        for i, score in enumerate(score_data):
            chromosome_dict.append([i, score])
            score_sum += score
        avg_score = score_sum/len(score_data)

        #get rid of genes that performed under the average and add multiples of
        # the good ones relative to how good they did
        final_selected = []
        for chrome in chromosome_dict:
            if chrome[1] >= avg_score:
                multiples = int(chrome[1]/avg_score)
                for i in range(multiples):
                    if len(final_selected) < self.genome_size:
                        final_selected.append(chrome)

        #fill up the remaining genome population with random chromosomes from
        # the previous genome.
        for i in range(self.genome_size - len(final_selected)):
            final_selected.append(chromosome_dict[random.randint(0, self.genome_size-1)])

        #rand = random.randint(0, int(score_sum))
        #random.shuffle(chromosome_dict)

        #return a list of the selected chromosomes
        return final_selected

    '''
    def single_sps_range(self, LEN):
        #2 section splitting for W2
        mr1=LEN-2
        mr2=LEN-1
        R1 = random.randint(0,LEN-1)
        R2 = random.randint(R1+1,LEN)
        return [R1,R2]

    def gen_rand_sps(self, LEN):
        #3 section splitting for W1
        mr1=LEN-2
        mr2=LEN-1
        mr3=LEN
        R1 = random.randint(0,mr1-1)
        R2 = random.randint(R1+1,mr1)
        R3 = random.randint(R2,mr2-1)
        R4 = random.randint(R3+1,mr2)
        R5 = random.randint(R4,mr3-1)
        R6 = random.randint(R5+1,mr3)
        return [[R1,R2],[R3,R4],[R5,R6]]
    '''
    def gene_swap(self, parA, parB):
        #Currently Universal Probability Operator Swapping - UPOS

        #W1_region_lst = self.gen_rand_sps(self.hidden_layer_size)
        #W2_region = self.single_sps_range(self.output_layer_size)
        parA_copy = parA.copy()
        parB_copy = parB.copy()

        #W1 gene traversal
        for bitstrip in range(self.input_size):
            for gene in range(self.hidden_layer_size):
                if self.prob_chance(self.crossover_rate):
                    a = parA["W1"][bitstrip][gene]
                    b = parB["W1"][bitstrip][gene]
                    parA_copy["W1"][bitstrip][gene] = b
                    parB_copy["W1"][bitstrip][gene] = a

        #W2 gene traversal
        for bitstrip in range(self.hidden_layer_size):
            for gene in range(self.output_layer_size):
                if self.prob_chance(self.crossover_rate):
                    a = parA["W2"][bitstrip][gene]
                    b = parB["W2"][bitstrip][gene]
                    parA_copy["W2"][bitstrip][gene] = b
                    parB_copy["W2"][bitstrip][gene] = a
        return [parA, parB]

    def crossover(self, selected_chromosomes):
        crossovers = self.selection_num
        current_species = self.get_current_species()
        current_genome = self.get_current_genome()
        genome = self.readFile(path.join(self.genetics_dir, "s_"+str(current_species)+"_genome_"+str(current_genome)+".json"))

        #Gathering selected chromosome genetic data and add them to the next genome
        #parent1 = genome["chromosome_"+str(selected_chromosomes[0][0])]
        mate_parents = []
        offspring = []
        #for i in range(self.selection_num):
        for i in range(len(selected_chromosomes)):
            mate_parents.append(genome["chromosome_"+str(selected_chromosomes[i][0])])
            offspring.append(genome["chromosome_"+str(selected_chromosomes[i][0])])
        print("TESTESTESTESTEST::::: ", mate_parents[0])
        #mate the selected parents and create 2 children for each mating and add them to the rest of the genome
        for i in range(0, len(mate_parents)):
            #parA = mate_parents[i-1]
            '''FLAG'''
            parA = mate_parents[random.randint(0, len(mate_parents)-1)]
            #print("len ", len(mate_parents))
            #print("i ", i)
            parB = mate_parents[i]
            children = self.gene_swap(parA, parB)
            offspring.extend(children)

        #cloning to fill the rest of the genome
        if self.cloning_enabled:
            for i in range(self.genome_size - len(offspring)):
                offspring.append(mate_parents[0])

        return offspring

    def prob_chance(self, rate):
        chance_pool = []
        for yes in range(int(rate)):
            chance_pool.append(1)
        for no in range(100-len(chance_pool)):
            chance_pool.append(0)
        random.shuffle(chance_pool)
        answer = chance_pool[random.randint(0, 99)]
        if answer == 1:
            return True
        elif answer == 0:
            return False

    def mutation(self, offspring_chromosomes):
        '''
        take the offspring and "mutate" them by randomly selecting a random
        size selection in a chromosome and swap the gene values to their
        opposite values. 0 -> 1 & 1 -> 0

        Parameters: ( (list of dict) offspring_chromosomes ) should be size 6
        '''
        #TODO: make this method a lot faster. Might have to merge with crossover method.
        mutated_offspring = {}
        #offspring_chromosomes is a list of all of the chromosomes
        for i, offspring in enumerate(offspring_chromosomes):
            off_chrome = offspring.copy()
            #W1 gene traversal
            for bitstrip in range(self.input_size):
                for gene in range(self.hidden_layer_size):
                    if self.prob_chance(self.mutation_rate):
                        if self.activation_func == "SIGMOID":
                            #for sigmoid
                            off_chrome["W1"][bitstrip][gene] = 1 - off_chrome["W1"][bitstrip][gene]
                        elif self.activation_func == "TANH":
                            #For use with tanh activation function
                            off_chrome["W1"][bitstrip][gene] = -1 * off_chrome["W1"][bitstrip][gene]

            #W2 gene traversal
            for bitstrip in range(self.hidden_layer_size):
                for gene in range(self.output_layer_size):
                    if self.prob_chance(self.mutation_rate):
                        if self.activation_func == "SIGMOID":
                            #for sigmoid
                            off_chrome["W2"][bitstrip][gene] = 1 - off_chrome["W2"][bitstrip][gene]
                        elif self.activation_func == "TANH":
                            #For use with tanh activation function
                            off_chrome["W2"][bitstrip][gene] = -1 * off_chrome["W2"][bitstrip][gene]
            mutated_offspring[i] = off_chrome
        return mutated_offspring

    def update_reference_file(self, reference_key, summand, flag=0):
        #check if the reference file even exists, if not then initialize it.
        if not self.initialize_ref_file():
            rf = self.readFile(self.reference_file)
            if flag == 0:
                rf[reference_key] += summand
            elif flag == 1:
                rf[reference_key] = summand
            self.writeFile(self.reference_file, rf)
            #logging.info("Successfully updated reference file.")
        #else:
            #logging.warning("reference file was empty when trying to update. Initialized.")

    def update_genome_species(self, filename):
        #check if the reference file even exists, if not then initialize it.
        if not file_is_empty(filename):
            rf = self.readFile(filename)
            rf["species"] += 1
            self.writeFile(filename, rf)
        #else:
            #logging.warning("Could not update species. File is empty.")

    def new_score_file(self, data):
        #check if file exists
        species_idx = self.get_current_species()
        #the space after the last genome in the species
        genome_idx  = self.get_current_genome()
        newGenomeScoreFile = path.join(self.score_dir, "s_"+str(species_idx)+"_genome_"+str(genome_idx)+".json")
        try:
            os.stat(self.score_dir)
            #logging.info("Score Directory already created")
        except:
            os.mkdir(self.score_dir)
            #logging.info("No Score Directory. Created it.")

        if not path.exists(newGenomeScoreFile):
            new_data = {
                "species" : species_idx,
                "scores" : data
            }
            self.createFile(newGenomeScoreFile)
            self.writeFile(newGenomeScoreFile, new_data)
        #else:
            #then the file has already been created and there is a problem
            #logging.warning("Attempted to create new Genome Score file, but the given score file has already been created.")

    def new_genome_file(self, species_idx, genome_idx, data):
        newGenomeGeneFile = path.join(self.genetics_dir, "s_"+str(species_idx)+"_genome_"+str(genome_idx)+".json")
        try:
            os.stat(self.genetics_dir)
            #logging.info("genetics Directory already created")
        except:
            os.mkdir(self.genetics_dir)
            #logging.info("No Score Directory. Created it.")

        if not path.exists(newGenomeGeneFile):
            new_data = {
                "species" : species_idx,
            }
            for chromosome in data.keys():
                new_data[chromosome] = data[chromosome]

            #update reference file and add 1 to genome count
            self.update_reference_file("num_of_genomes", 1)
            self.update_reference_file("last_genome", 1)

            self.createFile(newGenomeGeneFile)
            self.writeFile(newGenomeGeneFile, new_data)
        #else:
            #then the file has already been created and there is a problem
            #logging.warning("Attempted to create new Genome Gene file, but given score file has already been created.")

    def get_current_species(self):
        '''
        Returns the last species added to the population.
        '''
        refData = self.readFile(self.reference_file)
        current_species = refData["num_of_species"]
        return current_species

    def get_current_genome(self):
        '''
        Returns the last genome added to the current_species.
        '''
        refData = self.readFile(self.reference_file)
        current_genome = refData["last_genome"]
        return current_genome

    def get_current_chromosome(self, flag=0, s_idx=0, g_idx=0, last_chromosome=0):
        '''
        Returns the gene data of a chromosome.

        Parameters: (int chromosome_idx)
        '''
        refData = self.readFile(self.reference_file)

        if flag == 1:
            return last_chromosome
        elif flag != 2:
            s_idx = self.get_current_species()
            g_idx = self.get_current_genome()
            last_chromosome = refData["last_chromosome"]
        genetics_file = self.readFile(path.join(self.genetics_dir, "s_"+str(s_idx)+"_genome_"+str(g_idx)+".json"))

        return genetics_file["chromosome_"+str(last_chromosome)]

    def new_chromosome(self):
        W1, W2 = 0,0
        if self.activation_func == "SIGMOID":
            #for sigmoid
            W1 = np.random.sample((self.input_size, self.hidden_layer_size))
            W2 = np.random.sample((self.hidden_layer_size, self.output_layer_size))
        elif self.activation_func == "TANH":
            #For use with tanh activation function
            W1 = 2 * np.random.sample((self.input_size, self.hidden_layer_size)) - 1
            W2 = 2 * np.random.sample((self.hidden_layer_size, self.output_layer_size)) - 1

        return W1, W2

    def new_genome(self, data = None):
        '''
        Takes a dictionary of offspring chromosomes fills the rest of them with random chromosomes
        then writes it to the genetics.json file. If no offspring data is given then it generates
        a dictionary of random chromosomes and writes that to the genetics.json file instead.

        Parameters: (dict data)
        '''
        to_species_idx = self.get_current_species()
        #the space after the last genome in the species
        to_genome_idx  = self.get_current_genome() + 1
        #number of chromosomes per genome
        size = self.genome_size

        new_data = {}
        if data == None:
            #add size amount of chromosomes
            for i in range(size):
                W1, W2 = self.new_chromosome()
                new_data["chromosome_"+str(i)] = { "W1" : W1.tolist(), "W2" : W2.tolist() }
        else:
            #set new_data equal to the offspring
            for i, chromosome_key in enumerate(data.keys()):
                new_data["chromosome_"+str(i)] = data[chromosome_key]

            #add random chromosomes in the remaining slots
            for i in range(size-len(data)):
                W1, W2 = self.new_chromosome()
                new_data["chromosome_"+str(len(data)+i)] = { "W1" : W1.tolist(), "W2" : W2.tolist() }

        #write new genome containing the series of chromosomes to json file
        self.update_reference_file("num_of_chromosomes", size)
        self.new_genome_file(to_species_idx, to_genome_idx, new_data)

    def new_species(self):
        '''
        checks the last two genomes for speciation and moves the last genome to a new species in
        both genetics.json and scores.json files
        '''
        current_species_idx = self.get_current_species()
        current_genome_idx = self.get_current_genome()

        if self.speciation_check():
            #speciation has occured
            #increment the species number in the reference file
            self.update_reference_file("num_of_species", 1)
            self.update_reference_file("last_genome", 0, flag=1)
            #change the "species" value in the current genome file
            current_score_file = path.join(self.score_dir, "s_"+str(current_species_idx)+"_genome_"+str(current_genome_idx)+".json")
            current_genetics_file = path.join(self.genetics_dir, "s_"+str(current_species_idx)+"_genome_"+str(current_genome_idx)+".json")
            self.update_genome_species(current_score_file)
            self.update_genome_species(current_genetics_file)

            #change the species number of the current genome and rename the file
            os.rename(current_score_file,
            path.join(self.score_dir, "s_"+str(current_species_idx + 1)+"_genome_"+str(0)+".json"))

            os.rename(current_genetics_file,
            path.join(self.genetics_dir, "s_"+str(current_species_idx + 1)+"_genome_"+str(0)+".json"))
        #else:
            #speication has not occured
            #do nothing
            #logging.warning("Speciation has not occurred between s_{}_genome_{} and s_{}_genome_{}.".format(current_species_idx,current_genome_idx-1 ,current_species_idx, current_genome_idx))

    def speciation_check(self, flag=0):
        '''
        compares the standard deviation of two genomes and checks if the percent difference is
        greater than 50%, if it is then create a new species. Shouldnt be in this part of the code.
        '''

        current_species = self.get_current_species()
        B_old_idx = self.get_current_genome() - 1
        A_New_idx = self.get_current_genome()
        try:
            genome_B_Old = self.readFile(path.join(self.score_dir, "s_"+str(current_species)+"_genome_"+str(B_old_idx)+".json"))
            genome_A_New = self.readFile(path.join(self.score_dir, "s_"+str(current_species)+"_genome_"+str(A_New_idx)+".json"))
        except:
            #logging.warning("Checking speciation after already checked or not enough genomes to check.")
            return False

        gA_scores = genome_B_Old["scores"]
        gB_scores = genome_A_New["scores"]

        SA = standard_deviation(gA_scores)
        SB = standard_deviation(gB_scores)
        PD = percent_difference(SA, SB)
        if flag == 0:
            if (PD >= 20) and (SB > SA):
                return True
        elif flag == 1:
            return PD
        return False

    def pop_convergence(self, flag=0):
        '''
        take the standard deviation of each species and calculate each percent difference until it
        approaches 0.
        '''
        #TODO: Fix this method. throws errors when on. needs further testing.
        min_approach = 0

        ref_file_data = self.readFile(self.reference_file)
        num_of_species = ref_file_data["num_of_species"]
        num_of_genomes = ref_file_data["num_of_genomes"]
        #score_data = self.readFile(self.score_file)
        #species_keys = data["population"].keys()

        if num_of_species > 1:
            species_SD_list = []
            current_species_SDs = []
            '''
            for species_key in species_keys:
                for genome_key in score_data[species_key].keys():
                    #take the standard deviation of each genome and add it to the species list
                    genome = data["population"][species_key][genome_key]
                    '''
            #loop through each genome file and take the standard deviation of each genome socre list and add it to the species list
            #-->     current_species_SDs.append(standard_deviation(genome))
                #take the standard deviation of each species and add it to the total species Sd list
                #species_SD_list.append(standard_deviation(current_species_SDs))

            for species_idx in range(num_of_species):
                for genome_idx in range(num_of_genomes):
                    #loop through each genome file and take the standard deviation of each genome socre list and add it to the species list
                    try:
                        genome_scores = self.readFile(path.join(self.score_dir,
                        "s_"+str(species_idx)+"_genome_"+str(genome_idx)+".json"))["scores"]
                        current_species_SDs.append(standard_deviation(genome_scores))
                    except:
                        #logging.info("Checking convergence for genome thats not there. moving on to next genome or species.")
                        continue

                #take the standard deviation of each species and add it to the total species Sd list
                species_SD_list.append(standard_deviation(current_species_SDs))
                #current_species_SDs = []

            #find the percent difference between each species
            population_diffs = []
            #species_SD_list.pop(0)
            #print(species_SD_list)
            for i in range(1, len(species_SD_list)):
                population_diffs.append(percent_difference(species_SD_list[i], species_SD_list[i-1]))

            #find the minimum number out of the percent differences
            #and if it is less than min_approach return true
            min_PD = 0
            #print(species_SD_list)
            #print(current_species_SDs)
            #print(population_diffs)
            try:
                min_PD = min(population_diffs)
            except:
                min_PD = 0
        if flag == 0:
            if min_PD <= min_approach:
                return True
            return False
        elif flag == 1:
            return min_PD

    def initialize_ref_file(self):
        if not path.exists(self.reference_file):
            self.createFile(self.reference_file)
        #else:
            #logging.warning("reference.json already exists, no need to recreate it.")

        if file_is_empty(self.reference_file):
            #set up initial file values
            data = {
                "num_of_species" : 0,
                "num_of_genomes" : 0,
                "num_of_chromosomes" : 0,
                "last_chromosome" : 0,
                "last_genome" : -1
            }

            self.writeFile(self.reference_file, data)
            return True
        else:
            #logging.warning("reference.json file is not empty. Did not initialize.")
            return False

    def readFile(self, filename):
        #reads data from a json file and returns it
        with open(filename, 'r') as f:
            data = json.load(f)
            return data
        f.close()

    def writeFile(self,filename, data):
        #writes data to a json file
        with open(filename, 'w') as f:
            json.dump(data, f, sort_keys=True, indent=1)
        f.close()

    def createFile(self, filename):
        f = open(filename, "x")
        #json.dump(data, f, sort_keys=True, indent=1)
        f.close()

##############_TESTING_##########################################_TESTING_##########################
NN = Neural_network()
#NN.update_reference_file(0,0,flag=1)
#NN.new_genome()
#NN.new_score_file([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
import time
t0 = time.time()
#NN.update_reference_file("last_chromosome", 0, flag=1)
#NN.new_score_file([999,555,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
selected = NN.selection_elitism()
bread = NN.crossover(selected)
mutated = NN.mutation(bread)
#NN.new_genome(mutated)
tf = time.time()
print(tf-t0)

print("-----")
NN.writeFile("pooties.json", mutated)
