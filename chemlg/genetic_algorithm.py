from __future__ import print_function
from builtins import range

from deap import base, creator, tools
import random
import pandas as pd
import time
import math

from chemlg.libgen import molecule, check_building_blocks, get_index_list
import pybel



class building_blocks(object):
    """A class for each of the building blocks that are read from the building blocks file.
    Class Variables are created for smiles, molecule code, number of atoms, list of indices of potential sites for reaction and the number of such sites.

    Parameters
    ----------
    mol: dict
        dict object created from the molecule function

    Returns
    -------

    """
    def __init__(self, mol):
        self.smiles = mol['reverse_smiles']
        self.smiles_struct = mol['code']
        mol_ob = pybel.readstring("smi", mol['reverse_smiles'])
        self.atom_len = len(list(mol_ob.atoms))
        self.index_list = get_index_list(mol, list(mol_ob.atoms)) 
        self.spaces = len(self.index_list)


class GeneticAlgorithm(object):
    """
            A genetic algorithm class for search or optimization problems, built on top of the
            Distributed Evolutionary Algorithms in Python (DEAP) library. There are three algorithms with different genetic
            algorithm selection methods. The documentation for each method is mentioned in the documentation for the search module.

            Parameters
            ----------
            evaluate: function
                The objective function that has to be optimized. The first parameter of the objective function should be
                an object of type <chromosome_type>. Objective function should return a tuple.

            bb_file: list
                building blocks file

            max_indi_len: int, optional (default=20)
                Integer to specify the maximum number of building blocks (similar or dissimilar) that can be used to construct a molecule. 

            weights: tuple of integer(s), optional (default = (-1.0, ) )
                A tuple containing fitness objective(s) for objective function(s). Ex: (1.0,) for maximizing and (-1.0,)
                for minimizing a single objective function

            pop_size: integer, optional (default = 50)
                An integer which denotes the size of the population

            n_generations: integer, optional (default = 20)
                An integer for the number of generations for evolving the population

            conv_criteria: int, optional (default=10)
                Integer specifying the maximum number of generations for which the algorithm can select the same best individual, after which 
                the search terminates.

            algorithm: int, optional (default=1)
                The algorithm to use for the search. Algorithm descriptions are in the documentation for the search method.

            """

    def __init__(self, 
                evaluate, 
                bb_file,
                max_indi_len = 20,
                weights=(-1.0, ), 
                pop_size=50, 
                n_generations=20, 
                conv_criteria=10, 
                algorithm=1):

        self.Weights = weights
        self.bb_file = bb_file
        self.max_ind_len = max_indi_len
        self.evaluate = evaluate
        self.pop_size = pop_size
        self.n_generations = n_generations
        self.conv_criteria = conv_criteria
        self.algo = algorithm
        
    def chromosome_generator(self):
        """Generates the chromosome for the algorithm, after reading and validating the molecules from the building blocks file. 

        Parameters
        ----------

        Returns
        -------

        """   
        ## Reading the building blocks from the input file
        initial_mols = []
        try :
            infile = open(self.bb_file)
        except:
            tmp_str = "Building blocks file "+self.bb_file+" does not exist. "
            tmp_str = tmp_str+"Please provide correct building blocks file.\n"
            exit(code=1)
            
        for i,line in enumerate(infile):
            smiles = line.strip()
            if smiles.isspace() or len(smiles)==0 or smiles[0]=='#':
                continue
            if '[X]' in smiles:
                smiles = smiles.replace('[X]','[Ra]')
            smiles = check_building_blocks(smiles,i+1,self.bb_file)

            # removing duplicates in the input list based on canonical smiles
            temp = molecule(smiles, 'F'+str(len(initial_mols)+1))
            is_duplicate = False
            for z in initial_mols:
                if temp['can_smiles'] not in z['can_smiles']:
                    continue
                is_duplicate = True
            if not is_duplicate:
                initial_mols.append(temp)
                
        # create building blocks class objects for each validated molecule and store them in a list. 
        self.bb = [building_blocks(i) for i in initial_mols]
        
        i = 0
        chromosome = []
        ind_len = random.randint(2, self.max_ind_len)
        while i < ind_len:
            if i == 0:
                r = random.randint(0, len(self.bb) - 1)             # randomly select building block
                chromosome.append(self.bb[r].smiles_struct)
                for j in range(self.bb[r].spaces):
                    chromosome.append([])
            else:
                avl_pos = count_list(chromosome)[0]
                if len(avl_pos) <= 0:
                    return chromosome
                r = random.randint(0, len(avl_pos) - 1)                 # random number for selecting handle of 1st bb
                s = random.randint(0, len(self.bb) - 1)             # random number for selecting bb
                t = random.randint(1, self.bb[s].spaces)            # random number for selecting handle of 2nd bb
                nested_lookup(chromosome, avl_pos[r]).append(self.bb[s].smiles_struct)
                for j in range(self.bb[s].spaces):
                    if (j+1) != t:
                        nested_lookup(chromosome, avl_pos[r]).append([])
                    else:
                        nested_lookup(chromosome, avl_pos[r]).append(['C'])
                
            i += 1

        return chromosome

    def pre_eval(self, indi_list):
        """Pre-processes the individuals/chromosomes before sending them to the evaluate function. The function converts the lists of lists generated by 
        the algorithm to actual molecules and then sends the molecules to the evaluate function.

        Parameters
        ----------
        indi_list: list,
            individual received from the algorithm.

        Returns
        -------
        fit_val: float,
            fitness value of the individual

        """
        indi_list = tuple(indi_list)
        mol_code_list, parent_list, handle_list, mol_combi = [], [], [], ''
        f_lists = count_list(indi_list)[1]

        parent_list.append([indi_list, indi_list, -1])
        while len(parent_list)!= 0:
            iterate_over = parent_list[-1][1]
            while len(iterate_over) != 0:
                for k_ind, k in enumerate(iterate_over):
                    # continue for loop if item already traversed
                    if k_ind <= parent_list[-1][2]:
                        continue

                    # track molecule smiles, its length and the list's nested index in the given list
                    if not isinstance(k, list):
                        for mol_objs in self.bb:
                            if mol_objs.smiles_struct == k:
                                new_mol_smi = mol_objs.smiles
                                mol_len = mol_objs.atom_len
                                break
                        mol_combi = mol_combi + new_mol_smi + '.'
                        for fl in f_lists:
                            if nested_lookup(indi_list, fl) == iterate_over:
                                l_ind = fl
                                break
                        if iterate_over == indi_list: l_ind = -1
                        mol_code_list.append((k, mol_len, l_ind))
                        
                    else:
                        if k:
                            if k[0] == 'C':
                                handle_2 = k_ind
                                handle_1 = parent_list[-2][2]
                                for mol_objs in self.bb:
                                    if mol_objs.smiles_struct == parent_list[-1][0][0]:
                                        handle_1 = mol_objs.index_list[handle_1]
                                    if mol_objs.smiles_struct == parent_list[-1][1][0]:
                                        handle_2 = mol_objs.index_list[handle_2]
                                ci_size = 0
                                # get the index numbers for both handles by checking the molecule codes and their list indices.
                                for mcl in mol_code_list:
                                    if mcl[0] == parent_list[-1][0][0] and nested_lookup(indi_list, mcl[2]) == parent_list[-1][0]:
                                        handle_1 += ci_size
                                    if mcl[0] == parent_list[-1][1][0] and nested_lookup(indi_list, mcl[2]) == parent_list[-1][1]:
                                        handle_2 += ci_size
                                        break
                                    ci_size += mcl[1]
                                # append the handle indices to a list
                                handle_list.append([handle_1, handle_2])
                            else:
                                parent_list[-1][2] = k_ind
                                parent_list.append([iterate_over, k, -1])
                                break
                if len(iterate_over) != 0:
                    del parent_list[-1]
                    break
        # read the collected smiles into pybel
        mol_combi = pybel.readstring('smi', mol_combi)
        # create bonds for each of the handles in handle_list
        for handles in handle_list:
            mol_combi.OBMol.AddBond(handles[0], handles[1], 1)
        mol_combi = mol_combi.write("can")
        mol_combi = pybel.readstring("smi", mol_combi)
        # todo: check for validity against generation rules.
        # if not if_add: fit_val = (-1000.0 * self.Weights[0],)


        # fit_val should be a tuple!! 
        fit_val = self.evaluate(mol_combi)
        return fit_val,

    def fit(self):
        """
        Setting up the DEAP - genetic algorithm parameters and functions.

        Notes
        -----
        Must call this function before calling the search method.
        """
        
        creator.create("FitnessMin", base.Fitness, weights=self.Weights)
        creator.create("Individual", list, fitness=creator.FitnessMin)
        self.toolbox = base.Toolbox()
        self.toolbox.register("individual", tools.initIterate,
                              creator.Individual, self.chromosome_generator)
        self.toolbox.register("population", tools.initRepeat, list,
                              self.toolbox.individual)

        self.toolbox.register("selectRoulette", tools.selRoulette)        
        self.toolbox.register("evaluate", self.pre_eval)
        
    def crossover(self, child1, child2):

        c1 = count_list(child1)[1]
        c2 = count_list(child2)[1]
        if not (len(c1)==0 or len(c2)==0):
            r1 = random.randint(0, len(c1) - 1)
            r2 = random.randint(0, len(c2) - 1)
            t1 = nested_lookup(child1, c1[r1])             # pick and store a list randomly from child 1
            t2 = nested_lookup(child2, c2[r2])             # pick and store a list randomly from child 2
            if len(c1[r1]) > 1:
                del nested_lookup(child1, c1[r1][:-1])[c1[r1][-1]]
                nested_lookup(child1, c1[r1][:-1]).insert(c1[r1][-1], t2)
            else:
                del child1[c1[r1][0]]
                child1.insert(c1[r1][0], t2)

            if len(c2[r2]) > 1:
                del nested_lookup(child2, c2[r2][:-1])[c2[r2][-1]]
                nested_lookup(child2, c2[r2][:-1]).insert(c2[r2][-1], t1)
            else:
                del child2[c2[r2][0]]
                child2.insert(c2[r2][0], t1)

    def custom_mutate(self, indi):
        t = random.random()
        if t < 0.33:                                                            # add a block randomly
            c = count_list(indi)[0]
            if not c:
                self.custom_mutate(indi)
            else:
                r = random.randint(0, len(c) - 1)                                   # random number for selecting empty handle
                s = random.randint(0, len(self.bb)-1)                           # random number for selecting which bb to insert
                t = random.randint(1, self.bb[s].spaces)                        # random number for selecting which handle to connect to in the new bb
                nested_lookup(indi, c[r]).append(self.bb[s].smiles_struct)
                for j in range(self.bb[s].spaces):
                    if (j+1) != t:
                        nested_lookup(indi, c[r]).append([])
                    else:
                        nested_lookup(indi, c[r]).append(['C'])
        elif 0.33 < t < 0.66:                                                   # delete a block randomly
            c = count_list(indi)[1]
            r = random.randint(0, len(c) - 1)
            if len(c[r]) > 1:
                del nested_lookup(indi, c[r])[:]
            else:
                self.custom_mutate(indi)
        else:                                                                   # replace a block randomly
            c = count_list(indi)[1]
            r = random.randint(0, len(c) - 1)
            while isinstance(nested_lookup(indi, c[r])[0], list):
                c[r].append(0)
            old_block_code = nested_lookup(indi, c[r])[0]
            for x in self.bb:
                if x.smiles_struct == old_block_code: old_block = x
            s = random.randint(0, len(self.bb)-1)
            new_block = self.bb[s]
            nested_lookup(indi, c[r])[0] = new_block.smiles_struct
            diff = old_block.spaces - new_block.spaces
            if diff < 0:
                for k in range(diff):
                    nested_lookup(indi, c[r]).append([])
            elif diff > 0:
                del nested_lookup(indi, c[r])[-diff:]

    def search(self, init_pop_frac = 0.35, crossover_pop_frac = 0.35):
        """
        Algorithm 1:
            Initial population is instantiated. 
            Roulette wheel selection is used for selecting individuals for crossover and mutation.
            The initial population, crossovered and mutated individuals form the pool of individuals from which the best
            n members are selected as the initial population for the next generation, where n is the size of population.

        Algorithm 2:
            Initial population is instantiated.
            Roulette wheel selection is used for selecting individuals for crossover and mutation.
            The initial population, crossovered and mutated individuals form 3 different pools of individuals. Based on
            input parameters 1 and 2, members are selected from each of these pools to form the initial population for the
            next generation. Fraction of mutated members to select for next generation is decided based on the two input
            parameters and the size of initial population.

        Algorithm 3:
            Initial population is instantiated.
            Roulette wheel selection is used for selecting individuals for crossover and mutation.
            The initial population, crossovered and mutated individuals form the pool of individuals from which n members
            are selected using Roulette wheel selection, but without replacement to ensure uniqueness of members in the next
            generation, as the initial population for the next generation, where n is the size of population.


        Parameters
        ----------
        init_pop_frac: float, optional (default = 0.4)
            Fraction of initial population to select for next generation

        crossover_pop_frac: float, optional (default = 0.3)
            Fraction of crossover population to select for next generation

        
        Returns
        -------
        best_ind_df:  pandas dataframe
            A pandas dataframe of best individuals of each generation

        best_ind:  list of <chromosome_length> numbers
            The best individual after the last generation.
        Notes
        -----
        """

        pop = self.toolbox.population(n=self.pop_size)
        # Evaluate the initial population
        fitnesses = map(self.toolbox.evaluate, pop)
        for ind, fit in zip(pop, fitnesses):
            ind.fitness.values = fit

        # Generate and evaluate crossover and mutation population
        if self.algo == 3:
            co_pop = []
            while len(co_pop) < int(math.ceil(0.8*len(pop))):
                c = self.toolbox.selectRoulette(pop, 1)
                if c not in co_pop:
                    co_pop = co_pop + c
        else:
            co_pop = self.toolbox.selectRoulette(pop, int(math.ceil(0.8*len(pop))))
        co_pop = list(map(self.toolbox.clone, co_pop))
        mu_pop = self.toolbox.selectRoulette(pop, int(math.ceil(0.3*len(pop))))
        mu_pop = list(map(self.toolbox.clone, mu_pop))

        for child1, child2 in zip(co_pop[::2], co_pop[1::2]):
            self.crossover(child1, child2)
            del child1.fitness.values
            del child2.fitness.values

        for mutant in mu_pop:
            self.custom_mutate(mutant)
            del mutant.fitness.values

        # Evaluate the crossover and mutated population
        if self.algo == 2:
            a, b = int(math.ceil(init_pop_frac*len(pop))), int(math.ceil(crossover_pop_frac*len(pop)))
            total_pop = tools.selBest(pop, a) + tools.selBest(co_pop, b) + tools.selBest(mu_pop, len(pop)-a-b)
        else:
            total_pop = pop + co_pop + mu_pop
        invalid_ind = [ind for ind in total_pop if not ind.fitness.valid]
        fitnesses = list(map(self.toolbox.evaluate, invalid_ind))
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        # fits = [indi.fitness.values[0] for indi in pop]

        best_indi_per_gen = []
        best_indi_fitness_values = []
        timer = []
        convergence = 0
        for g in range(self.n_generations):
            if convergence >= self.conv_criteria:
                print("The search converged with convergence criteria = ", self.conv_criteria)
                break
            else:
                st_time = time.time()
                # Select the next generation individuals
                offspring = tools.selBest(total_pop, self.pop_size)
                # Clone the selected individuals
                offspring = list(map(self.toolbox.clone, offspring))
                if self.algo == 3:
                    co_pop = []
                    while len(co_pop) < int(math.ceil(0.8*len(pop))):
                        c = self.toolbox.selectRoulette(pop, 1)
                        if c not in co_pop:
                            co_pop = co_pop + c
                else:
                    co_pop = self.toolbox.selectRoulette(offspring, int(math.ceil(0.8*len(pop))))
                co_pop = list(map(self.toolbox.clone, co_pop))
                mu_pop = self.toolbox.selectRoulette(offspring, int(math.ceil(0.3*len(pop))))
                mu_pop = list(map(self.toolbox.clone, mu_pop))

                for child1, child2 in zip(co_pop[::2], co_pop[1::2]):
                    self.crossover(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

                for mutant in mu_pop:
                    self.custom_mutate(mutant)
                    del mutant.fitness.values

                # Evaluate the crossover and mutated population
                if self.algo == 2:
                    a, b = int(math.ceil(init_pop_frac*len(pop))), int(math.ceil(crossover_pop_frac*len(pop)))
                    total_pop = tools.selBest(pop, a) + tools.selBest(co_pop, b) + tools.selBest(mu_pop, len(pop)-a-b)
                else:
                    total_pop = offspring + co_pop + mu_pop
                invalid_ind = [ind for ind in total_pop if not ind.fitness.valid]
                fitnesses = list(map(self.toolbox.evaluate, invalid_ind))
                for ind, fit in zip(invalid_ind, fitnesses):
                    ind.fitness.values = fit


                # Storing the best individuals after each generation
                best_individual = tools.selBest(total_pop, 1)[0]
                if len(best_indi_per_gen)>0 and best_individual==best_indi_per_gen[-1]:
                    convergence += 1
                best_indi_per_gen.append(list(best_individual))
                best_indi_fitness_values.append(best_individual.fitness.values[0])

                tot_time = (time.time() - st_time)/(60*60)
                timer.append(tot_time)

                b1 = pd.Series(best_indi_per_gen, name='Best_individual_per_gen')
                b2 = pd.Series(best_indi_fitness_values, name='Fitness_values')
                b3 = pd.Series(timer, name='Time')
                best_ind_df = pd.concat([b1, b2, b3], axis=1)

        print("\n \n Best Individuals of each generation are:  \n \n" , best_ind_df)
        print("\n \n Best individual after %s evolutions is %s " % (self.n_generations, best_individual))
        del creator.FitnessMin
        del creator.Individual
        return best_ind_df, best_individual


def count_list(l):
    """ Function to get the nested indices of empty and filled lists generated by genetic algorithm class.

    Parameters
    ----------
    l: list,
        nested list of list

    Returns
    -------
    e_list_index: list,
        the nested indices of empty lists
    f_list_index: list,
        the nested indices of filled/non-empty lists

    """
    
    e_list_index, f_list_index, iterate_over = [], [], []
    for e, g in zip(l, range(len(l))):
        if isinstance(e, list):
            if not e:
                temp = [g]
                e_list_index.append(temp)
            else:
                if e != ['C']:
                    temp = [g]
                    f_list_index.append(temp)
                    iterate_over.append(e)

    while len(iterate_over) != 0:
        f_list = []
        for x, prefactor in zip(iterate_over, f_list_index[-len(iterate_over):]):
            if not f_list_index:
                prefactor = []
            for e, g in zip(x, range(len(x))):
                if isinstance(e, list):
                    if e == x:                                  # this is to check for self-referential lists
                        e_list_index.append(prefactor + temp)
                    else:
                        if not e:
                            temp = [g]
                            e_list_index.append(prefactor + temp)
                        else:
                            if e != ['C']:
                                temp = [g]
                                f_list_index.append(prefactor + temp)
                                f_list.append(e)

        del iterate_over[:]
        for items in f_list:
            iterate_over.append(items)
    return e_list_index, f_list_index

def nested_lookup(n, idexs):
    """Function to fetch a nested sublist given its nested indices.

    Parameters
    ----------
    n: list,
        the main list in which to look for the sublist

    idexs: list,
        the indices of the sublist 

    Returns
    -------
    list: sublist with given indices

    """

    if len(idexs) == 1:
        return n[idexs[0]]
    return nested_lookup(n[idexs[0]], idexs[1:])

