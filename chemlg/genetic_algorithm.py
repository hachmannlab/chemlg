from __future__ import print_function
from builtins import range
import itertools
import random
import pandas as pd
import time
import math
from copy import deepcopy
from chemlg.libgen import molecule, check_building_blocks, get_index_list, if_add, get_rules
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

            bb_file: str,
                Path to the building blocks file

            config_file: str,
                Path to the config file for generating the library

            output_dir: str, default = './'
                Path to the output directory.

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
                config_file,
                output_dir = './',
                max_indi_len = 20,
                fitness="Max", 
                crossover_size=0.5,
                mutation_size=0.3,
                pop_size=50, 
                algorithm=1,
                initial_population=None):

        self.bb_file = bb_file
        self.config_file = config_file
        self.output_dir = output_dir
        self.max_ind_len = max_indi_len
        self.evaluate = evaluate
        self.pop_size = pop_size
        self.algo = algorithm
        self.population, self.fitness_dict = None, {}
        self.crossover_size = crossover_size
        self.mutation_size = mutation_size
        self.initial_pop = initial_population
        if fitness.lower() == 'max': self.fit_val = 1
        else: self.fit_val = -1

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
            smiles = check_building_blocks(smiles,i+1,self.bb_file,self.output_dir)

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

        # Reading rules file
        rulesFile = open(self.config_file)
        self.rules_dict, args = get_rules(rulesFile, self.output_dir)
                
    def pop_generator(self, n):
        pop = []
        if self.initial_pop is not None:
            for i in self.initial_pop:
                pop.append(tuple(i))
        else:
            for _ in range(n):
                pop.append(self.chromosome_generator())
        return pop

    def chromosome_generator(self):
        """Generates the chromosome for the algorithm, after reading and validating the molecules from the building blocks file. 

        Parameters
        ----------

        Returns
        -------

        """   
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
                s = random.randint(0, len(self.bb) - 1)                 # random number for selecting bb
                t = random.randint(1, self.bb[s].spaces)                # random number for selecting handle of 2nd bb
                nested_lookup(chromosome, avl_pos[r]).append(self.bb[s].smiles_struct)
                for j in range(self.bb[s].spaces):
                    if (j+1) != t:
                        nested_lookup(chromosome, avl_pos[r]).append([])
                    else:
                        nested_lookup(chromosome, avl_pos[r]).append(['C'])
                
            i += 1

        return chromosome

    def list_to_smiles(self, indi_list):
        """The function converts the lists of lists generated by the algorithm to actual molecules.

        Parameters
        ----------
        indi_list: list,
            individual received from the algorithm.

        Returns
        -------
        mol_combi: str,
            canonical smiles of the molecule
        """   
        mol_code_list, parent_list, handle_list, mol_combi, mol_len = [], [], [], '', 0
        f_lists = count_list(indi_list)[1]

        parent_list.append([indi_list, indi_list, -1])
        while len(parent_list)!= 0:
            iterate_over = parent_list[-1][1]
            new_item = False
            for k_ind, k in enumerate(iterate_over):
                # continue for loop if item already traversed
                if k_ind <= parent_list[-1][2]:
                    continue

                # track molecule smiles, its length and the list's nested index in the given list
                if not isinstance(k, list):
                    for mol_objs in self.bb:
                        if mol_objs.smiles_struct == k:
                            new_mol_smi = mol_objs.smiles
                            mol_len += mol_objs.atom_len
                            break
                    mol_combi = mol_combi + new_mol_smi + '.'
                    if iterate_over == indi_list: nested_ind = -1
                    else:
                        iterate_over.insert(0, 'AB')
                        for fl in f_lists:
                            if nested_lookup(indi_list, fl) == iterate_over:
                                nested_ind = fl
                                break
                        del iterate_over[0]
                    mol_code_list.append((k, mol_len, nested_ind))
                    
                else:
                    if k:
                        if k[0] == 'C':
                            handle_2 = k_ind
                            handle_1 = parent_list[-2][2]
                            for mol_objs in self.bb:
                                if mol_objs.smiles_struct == parent_list[-1][0][0]:
                                    handle_1 = mol_objs.index_list[handle_1 - 1]
                                if mol_objs.smiles_struct == parent_list[-1][1][0]:
                                    handle_2 = mol_objs.index_list[handle_2 - 1]
                            
                            # get the index numbers for both handles by checking the molecule codes and their list indices.
                            for ind, mcl in enumerate(mol_code_list):
                                if mcl[2] == -1:
                                    continue
                                if mcl[0] == parent_list[-1][0][0]:
                                    parent_list[-1][0][0] += 'WW'
                                    if 'WW' in nested_lookup(indi_list, mcl[2])[0]: 
                                        handle_1 += mol_code_list[ind-1][1]
                                    parent_list[-1][0][0] = parent_list[-1][0][0][:-2]
                                
                                if mcl[0] == parent_list[-1][1][0]:
                                    parent_list[-1][1][0] += 'XX'
                                    if 'XX' in nested_lookup(indi_list, mcl[2])[0]: 
                                        handle_2 += mol_code_list[ind-1][1]
                                    parent_list[-1][1][0] = parent_list[-1][1][0][:-2]
                                    
                            # append the handle indices to a list
                            handle_list.append([handle_1, handle_2])
                        else:
                            parent_list[-1][2] = k_ind
                            parent_list.append([iterate_over, k, -1])
                            new_item = True
                            break
            if not new_item:
                del parent_list[-1]
                
        # read the collected smiles into pybel
        mol_combi = pybel.readstring('smi', mol_combi)
        # create bonds for each of the handles in handle_list
        for handles in handle_list:
            mol_combi.OBMol.AddBond(handles[0], handles[1], 1)
        mol_combi = mol_combi.write("can")
        return mol_combi

    def pre_eval(self, indi_list):
        """Pre-processes the individuals/chromosomes before sending them to the evaluate function. 

        Parameters
        ----------
        indi_list: list,
            individual received from the algorithm.

        Returns
        -------
        fit_val: float,
            fitness value of the individual

        """
        
        mol_combi_smiles = self.list_to_smiles(indi_list)
        mol_combi = pybel.readstring("smi", mol_combi_smiles)

        if not if_add(mol_combi_smiles, self.rules_dict, code='AA'): 
            fit_val = -1000.0 * self.fit_val
        else:
            fit_val = self.evaluate(mol_combi)
        return fit_val
        
    def crossover(self, child1, child2, fitness_dict):
        cn = True
        while cn:
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
            if tuple(child1) not in fitness_dict.keys() and tuple(child2) not in fitness_dict.keys(): cn = False
        return deepcopy(child1), deepcopy(child2)

    def custom_mutate(self, indi, fitness_dict):
        cn = True
        while cn:
            t = random.random()
            if t < 0.33:                                                            # add a block randomly
                c = count_list(indi)[0]
                if not c:
                    self.custom_mutate(indi, fitness_dict)
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
                    self.custom_mutate(indi, fitness_dict)
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
            if indi not in fitness_dict.keys(): cn = False
        return indi

    def RouletteWheelSelection(self, population, fit_dict, num):
        o_fits = [fit_dict[i] for i in population]
        fitnesses = [(((i-min(o_fits))/(max(o_fits)-min(o_fits))) + 1) for i in o_fits]
        if self.fit_val == -1:
            fitnesses = [fit**self.fit_val for fit in fitnesses]
        total_fitness = float(sum(fitnesses))
        rel_fitness = [f/total_fitness for f in fitnesses]
        # Generate probability intervals for each individual
        probs = [sum(rel_fitness[:i+1]) for i in range(len(rel_fitness))]
        # Draw new population
        new_population = []
        for _ in range(num):
            r = random.random()
            for i, individual in enumerate(population):
                if r <= probs[i]:
                    new_population.append(deepcopy(individual))
                    break
        return new_population

    def selectbest(self, pop, n, fitness_dict):
        best = []
        o_fits = [fitness_dict[i] for i in pop]
        fitnesses = [(((i-min(o_fits))/(max(o_fits)-min(o_fits))) + 1) for i in o_fits]
        if self.fit_val == -1:
            fitnesses = [fit**self.fit_val for fit in fitnesses]
        fits_sort = sorted(fitnesses, reverse=True)
        for i in range(min(n, len(pop))):
            best.append(deepcopy(pop[fitnesses.index(fits_sort[i])]))
        return best

    def search(self, n_generations=20, early_stopping=10, init_ratio = 0.35, crossover_ratio = 0.35):
        """
        Algorithm 1:
            Initial population is instantiated. 
            Roulette wheel selection is used for selecting individuals for crossover and mutation.
            The initial population, crossovered and mutated individuals form the pool of individuals from which the best
            n members are selected as the initial population for the next generation, where n is the size of population.

        Algorithm 2:
            Same as algorithm 1 but when selecting individuals for next generation, n members are selected using Roulette wheel selection.

        Algorithm 3:
            Same as algorithm 1 but when selecting individuals for next generation, best members from each of the three pools (initital population, crossover and mutation) are selected according to the input parameters in the search method.

        Algorithm 4:
            Same as algorithm 1 but mutation population is selected from the crossover population and not from the parents directly.


        Parameters
        ----------
        n_generations: integer, optional (default = 20)
                An integer for the number of generations for evolving the population

        early_stopping: int, optional (default=10)
                Integer specifying the maximum number of generations for which the algorithm can select the same best individual, after which 
                the search terminates.

        init_ratio: float, optional (default = 0.4)
            Fraction of initial population to select for next generation. Required only for algorithm 3.

        crossover_ratio: float, optional (default = 0.3)
            Fraction of crossover population to select for next generation. Required only for algorithm 3.

        
        Returns
        -------
        best_ind_df:  pandas dataframe
            A pandas dataframe of best individuals of each generation

        best_ind:  dict,
            The best individual after the last generation.

        """
        def fit_eval(invalid_ind, fitness_dict):
            if invalid_ind: 
                invalid_ind = [i for i in invalid_ind if i not in fitness_dict.keys()]
                fitnesses = list(map(self.pre_eval, invalid_ind))
                for ind, fit in zip(invalid_ind, fitnesses):
                    fitness_dict[tuple(ind)] = fit
            return fitness_dict


        if init_ratio >=1 or crossover_ratio >=1 or (init_ratio+crossover_ratio)>=1: raise Exception("Sum of parameters init_ratio and crossover_ratio should be in the range (0,1)")
        if self.population is not None:
            pop = self.population
            fitness_dict = self.fitness_dict
        else:
            pop = self.pop_generator(n=self.pop_size)       # list of tuples
            fitness_dict = {}
        
        # Evaluate the initial population
        fitness_dict = fit_eval(pop, fitness_dict)

        best_indi_per_gen, best_indi_fitness_values, timer, total_pop, convergence = [], [], [], [], 0
        for _ in range(n_generations):
            if convergence >= early_stopping:
                print("The search converged with convergence criteria = ", early_stopping)
                break
            else:
                st_time = time.time()
                cross_pop, mutant_pop, co_pop, psum = [], [], [], len(list(fitness_dict.items()))
                # Generate crossover population
                co_pop = self.RouletteWheelSelection(pop, fitness_dict, int(math.ceil(self.crossover_size*len(pop))))
                co_pop = list(itertools.combinations(co_pop, 2))
                combi = list(itertools.combinations(list(set(pop + total_pop)), 2))
                for child1, child2 in co_pop + combi:
                    if (len(list(fitness_dict.items())) - psum) >= int(math.ceil(self.crossover_size*len(pop))): break
                    c1, c2 = self.crossover(child1, child2, fitness_dict)
                    fitness_dict = fit_eval([c1, c2], fitness_dict)
                    cross_pop.extend([c1, c2])
                    
                # Generate mutation population
                if self.algo == 4:
                    mu_pop = self.RouletteWheelSelection(cross_pop, fitness_dict, int(math.ceil(self.pop_size*self.mutation_size)))
                else:
                    mu_pop = self.RouletteWheelSelection(pop, fitness_dict, int(math.ceil(self.mutation_size*len(pop))))
                
                for mutant in mu_pop:
                    mutant_pop.append(self.custom_mutate(mutant, fitness_dict))
                fitness_dict = fit_eval(mutant_pop, fitness_dict)
                
                # Select the next generation individuals
                total_pop = pop + cross_pop + mutant_pop
                if self.algo == 2:
                    pop = self.RouletteWheelSelection(total_pop, fitness_dict, self.pop_size)
                elif self.algo == 3:
                    p1 = self.selectbest(pop, int(init_ratio*self.pop_size), fitness_dict)
                    p2 = self.selectbest(cross_pop, int(crossover_ratio*self.pop_size), fitness_dict)
                    p3 = self.selectbest(mutant_pop, self.pop_size-len(p1)-len(p2), fitness_dict)
                    pop = p1 + p2 + p3
                else: pop = self.selectbest(total_pop, self.pop_size, fitness_dict)
                
                # Storing the best individuals after each generation
                best_individual = pop[0]
                if len(best_indi_per_gen)>0:
                    if best_individual==best_indi_per_gen[-1]: convergence += 1
                    else: convergence = 0
                best_indi_per_gen.append(best_individual)
                best_indi_fitness_values.append(fitness_dict[best_individual])
                tot_time = (time.time() - st_time)/(60*60)
                timer.append(tot_time)
                b1 = pd.Series(best_indi_per_gen, name='Best_individual')
                b2 = pd.Series(best_indi_fitness_values, name='Fitness_values')
                b3 = pd.Series(timer, name='Time (hours)')
                best_ind_df = pd.concat([b1, b2, b3], axis=1)
    

        # self.population = fitness_dict    # stores best individuals of last generation
        self.population = pop    # stores best individuals of last generation
        self.fitness_dict = fitness_dict
        # best_ind_dict = {}
        # for name, val in zip(self.var_names, best_individual):
            # best_ind_dict[name] = val
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

