import random, time, json

l1 = {"W1":[[2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,2],
            [0,2,0,2,2,0,0,2,0,0,2,0,2,2,0,2],
            [0,2,0,2,2,0,0,2,0,0,2,0,2,2,0,2],
            [0,2,0,2,2,0,0,2,0,0,2,0,2,2,0,2],
            [0,2,0,2,2,0,0,2,0,0,2,0,2,2,0,2],
            [0,2,0,2,2,0,0,2,0,0,2,0,2,2,0,2]],
      "W2":[[2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0],
            [2,0,2,0]]
      }

l2 = {"W1":[[1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,1],
            [1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,1],
            [1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,1],
            [1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,1],
            [1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,1],
            [1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,1]],
      "W2":[[1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0],
            [1,0,1,0]]
      }

input_size = 6
hidden_layer_size = 16
output_layer_size = 4
crossover_rate = 50

def prob_chance(rate):
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

def gene_swap(parA, parB):
	#Currently Universal Probability Operator Swapping - UPOS

	#W1_region_lst = self.gen_rand_sps(self.hidden_layer_size)
	#W2_region = self.single_sps_range(self.output_layer_size)
	parA_copy = parA.copy()
	parB_copy = parB.copy()

	#W1 gene traversal
	for bitstrip in range(input_size):
		for gene in range(hidden_layer_size):
			if prob_chance(crossover_rate):
				a = parA["W1"][bitstrip][gene]
				b = parB["W1"][bitstrip][gene]
				parA_copy["W1"][bitstrip][gene] = b
				parB_copy["W1"][bitstrip][gene] = a

	#W2 gene traversal
	for bitstrip in range(hidden_layer_size):
		for gene in range(output_layer_size):
			if prob_chance(crossover_rate):
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

def writeFile(filename, data):
	#writes data to a json file
	with open(filename, 'w') as f:
		json.dump(data, f, sort_keys=True, indent=1)
	f.close()

offspring = gene_swap(l1,l2)

writeFile("test1.json", offspring)
#writeFile("test2.json", offspring[1])
print(offspring[0])




