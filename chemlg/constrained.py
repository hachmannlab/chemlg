from __future__ import print_function
from builtins import range
import itertools
import random
import pandas as pd
import time
import math
from copy import deepcopy
import os
from chemlg.libgen import *


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
                The objective function that has to be optimized. The first argument of the function is an pybel object of a molecule. Objective function should return a tuple of desired target properties.

            fitness: tuple of tuple(s),
                A tuple of tuples for describing desired target properties. For each target property, tuple contains 2 values, the required optima and the cutoff value. For 'max' optima, the cutoff is the lower acceptable limit. For 'min' optima, the cutoff is the maximum allowed value for that property. Ex: (('max', 5.6), ('min', 20))

            bb_file: str,
                Path to the building blocks file

            config_file: str,
                Path to the config file for generating the library

            output_dir: str, default = './'
                Path to the output directory.

            crossover_size: int, optional (default = 50)
                size of crossover population

            mutation_size: integer, optional (default = 50)
                size of mutation population. Sum of crossover and mutation population is the total size of population generated in each generation.

            algorithm: int, optional (default=1)
                The algorithm to use for the search. Algorithm descriptions are in the documentation for the search method.

            """

    def __init__(self, 
                evaluate,
                fitness, 
                bb_file,
                config_file,
                output_dir = './',
                crossover_size=50,
                mutation_size=50,
                algorithm=1):

        self.bb_file = bb_file
        self.config_file = config_file
        self.output_dir = output_dir
        self.evaluate = evaluate
        self.algo = algorithm
        self.population, self.fit_list = None, ()
        self.crossover_size = int(crossover_size)
        self.mutation_size = int(mutation_size)
        self.fit_val = []
        self.pop_size = self.crossover_size + self.mutation_size
        for i in fitness:
            if i[1] == 0: print_lle("Cutoff values in the fitness cannot be zero.", self.output_dir)
            if i[0].lower() == 'max': self.fit_val.append((1, i[1]))
            else: self.fit_val.append((-1, i[1]))

        txt = "\n\n\n============================================================================================================"
        txt += "\n ChemLG - A Molecular and Materials Library Generator for the Enumeration and Exploration of Chemical Space"
        txt += "\n============================================================================================================\n\n\n"
        txt += "Program Version: 0.2 \t\t Release Date: Feb 20, 2019\n\n"
        txt += "(C) 2015-2018 Johannes Hachmann, Mohammad Atif Faiz Afzal \nUniversity at Buffalo - The State University of New York (UB)\n"
        txt += "Contact: hachmann@buffalo.edu, m27@buffalo.edu \nhttp://hachmannlab.cbe.buffalo.edu\n\n"
        txt += "With contributions by: \nJanhavi Abhay Dudwadkar (UB): Jupyter GUI\n\n"
        txt += "ChemLG is based upon work supported by the U.S. National Science Foundation under grant #OAC-1751161. \nIt was also supported by start-up funds provided by UB's School of Engineering and Applied Science and \nUB's Department of Chemical and Biological Engineering, the New York State Center of Excellence in Materials Informatics \nthrough seed grant #1140384-8-75163, and the U.S. Department of Energy under grant #DE-SC0017193."
        txt += "\n\nChemLG is copyright (C) 2015-2018 Johannes Hachmann and Mohammad Atif Faiz Afzal, all rights reserved. \nChemLG is distributed under 3-Clause BSD License (https://opensource.org/licenses/BSD-3-Clause). \n\n\n"
        print_ll(txt, self.output_dir)
        print_ll("===================================================================================================", self.output_dir)
        
        # Reading rules file
        try :
            rulesFile = open(self.config_file)
        except:
            tmp_str = "Config file does not exist. "
            tmp_str = tmp_str+"Please provide correct config file.\n"
            print_lle(tmp_str, self.output_dir,"Aborting due to wrong file.")
        print_ll("Reading generation rules", self.output_dir)
        print_ll("===================================================================================================", self.output_dir)
        self.rules_dict, args = get_rules(rulesFile, self.output_dir, {'rank': 0})
        combi_type, gen_len, outfile_type, max_fpf, lib_name = args
        gen_len, max_fpf = int(gen_len), int(max_fpf)
        if gen_len == 0:
            self.rules_dict['bb_final_lib'] = True

        
        ## Reading the building blocks from the input file
        initial_mols, i_smi_list = [], []
        print_ll("===================================================================================================", self.output_dir)
        print_ll("Reading building blocks from the file "+self.bb_file, self.output_dir)
        print_ll("===================================================================================================", self.output_dir)
        try :
            infile = open(self.bb_file)
        except:
            tmp_str = "Building blocks file "+self.bb_file+" does not exist. "
            tmp_str += "Please provide correct building blocks file.\n"
            print_lle(tmp_str, self.output_dir,"Aborting due to wrong file.")
        
        for i,line in enumerate(infile):
            smiles = line.strip()
            if smiles.isspace() or len(smiles)==0 or smiles[0]=='#':
                continue
            if '[X]' in smiles:
                smiles = smiles.replace('[X]','[Ra]')
            smiles = check_building_blocks(smiles,i+1,self.bb_file,self.output_dir, {'rank': 0})

            # removing duplicates in the input list based on canonical smiles
            temp = molecule(smiles, 'F'+str(len(initial_mols)+1))
            is_duplicate = False
            for z in initial_mols:
                if temp['can_smiles'] not in z['can_smiles']:
                    continue
                is_duplicate = True
            if not is_duplicate:
                initial_mols.append(temp)
                i_smi_list.append(temp['can_smiles'])
                
        print_ll('Number of buidling blocks provided = '+str(len(initial_mols))+'\n', self.output_dir)
        print_ll('unique SMILES: ', self.output_dir)
        print_ll(i_smi_list, self.output_dir)
        # create building blocks class objects for each validated molecule and store them in a list. 
        self.bb = [building_blocks(i) for i in initial_mols]

    def pop_generator(self, n):
        pop = []
        for _ in range(n):
            pop.append(tuple(self.chromosome_generator()))
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
        ind_len = random.randint(2, 5)
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

        return deepcopy(chromosome)

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

        # parent list: [list, sublist, index of sublist in list]
        parent_list.append([indi_list, indi_list, -100])
        while len(parent_list)!= 0:
            # iterate over new sublist
            iterate_over = parent_list[-1][1]
            new_item = False
            for k_ind, k in enumerate(iterate_over):
                # continue for loop if item already traversed
                if k_ind <= parent_list[-1][2]: continue

                if not isinstance(k, list):
                    for mol_objs in self.bb:
                        if mol_objs.smiles_struct == k:
                            new_mol_smi = mol_objs.smiles
                            mol_len += mol_objs.atom_len
                            break
                    mol_combi = mol_combi + new_mol_smi + '.'
                    if iterate_over == indi_list: nested_ind = -50
                    else:
                        iterate_over.insert(0, 'AB')
                        for fl in f_lists:
                            if nested_lookup(indi_list, fl) == iterate_over:
                                nested_ind = fl
                                break
                        del iterate_over[0]
                    # mol_code_list: [current sublist's molecule code, cumulative sum of atoms of all building blocks encountered so far, nested indices of current sublist]
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
                                if mcl[2] == -50:
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
                            parent_list.append([iterate_over, k, -100])
                            new_item = True
                            break
            if not new_item:
                del parent_list[-1]
                
        # read the collected smiles into pybel
        mol_combi = pybel.readstring('smi', mol_combi)
        # create bonds for each of the handles in handle_list
        for handles in handle_list:
            mol_combi.OBMol.AddBond(handles[0], handles[1], 1)
        for atom in list(mol_combi.atoms):
            ## Removing Francium and Radium atoms
            if atom.OBAtom.GetAtomicNum() == 87 or atom.OBAtom.GetAtomicNum() == 88:
                atom.OBAtom.SetAtomicNum(1)
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
        mol_combi_smiles = self.list_to_smiles(deepcopy(list(indi_list)))
        mol_combi = pybel.readstring("smi", mol_combi_smiles)

        if not if_add(mol_combi_smiles, self.rules_dict, code='a'): 
            return mol_combi_smiles, None
        else:
            fit_val = self.evaluate(mol_combi)
            if isinstance(fit_val, (tuple, list)): return mol_combi_smiles, tuple(fit_val)
            else: return mol_combi_smiles, tuple([fit_val])
        
    def crossover(self, child1, child2):
        child1, child2 = deepcopy(list(child1)), deepcopy(list(child2))
        c1 = count_list(child1)[1]
        c2 = count_list(child2)[1]
        if not (len(c1)==0 or len(c2)==0):
            r1 = random.randint(0, len(c1) - 1)
            r2 = random.randint(0, len(c2) - 1)
            t1 = nested_lookup(child1, c1[r1])                          # pick and store a list randomly from child 1
            t2 = nested_lookup(child2, c2[r2])                          # pick and store a list randomly from child 2
            if len(c1[r1]) > 1:
                del nested_lookup(child1, c1[r1][:-1])[c1[r1][-1]]
                nested_lookup(child1, c1[r1][:-1]).insert(c1[r1][-1], deepcopy(t2))
            else:
                del child1[c1[r1][0]]
                child1.insert(c1[r1][0], deepcopy(t2))

            if len(c2[r2]) > 1:
                del nested_lookup(child2, c2[r2][:-1])[c2[r2][-1]]
                nested_lookup(child2, c2[r2][:-1]).insert(c2[r2][-1], deepcopy(t1))
            else:
                del child2[c2[r2][0]]
                child2.insert(c2[r2][0], deepcopy(t1))
        return tuple(deepcopy(child1)), tuple(deepcopy(child2))

    def custom_mutate(self, indi):
        indi = deepcopy(list(indi))
        t = ['add', 'del', 'replace']
        random.shuffle(t)
        for i in t:
            if i == 'add':                                                          # add a block randomly
                c = count_list(indi)[0]
                if not c:
                    continue
                else:
                    r = random.randint(0, len(c) - 1)                               # random number for selecting empty handle
                    s = random.randint(0, len(self.bb)-1)                           # random number for selecting which bb to insert
                    t = random.randint(1, self.bb[s].spaces)                        # random number for selecting which handle to connect to in the new bb
                    nested_lookup(indi, c[r]).append(self.bb[s].smiles_struct)
                    for j in range(self.bb[s].spaces):
                        if (j+1) != t:
                            nested_lookup(indi, c[r]).append([])
                        else:
                            nested_lookup(indi, c[r]).append(['C'])
                    break
            elif i == 'del':                                                        # delete a block randomly
                c = count_list(indi)[1]
                r = random.randint(0, len(c) - 1)
                if not c or len(c[r]) < 2: continue
                del nested_lookup(indi, c[r])[:]
                break
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
                    for _ in range(-diff):
                        nested_lookup(indi, c[r]).append([])
                elif diff > 0:
                    tmp = deepcopy(nested_lookup(indi, c[r])[1:])
                    del nested_lookup(indi, c[r])[1:]
                    nested_lookup(indi, c[r]).append(['C'])
                    for i in range(new_block.spaces-1): nested_lookup(indi, c[r]).append(random.choice(tmp))
                                                    
                break
        return tuple(deepcopy(indi))

    def select(self, population, fit_list, num, choice="Roulette"):
        epop, efits = [i[0] for i in fit_list], [i[1] for i in fit_list]
        o_fits = [efits[epop.index(i)] for i in population]

        df_fits = pd.DataFrame(o_fits)
        # calculate distance from cut-offs
        weights = pd.DataFrame([df_fits[i] / self.fit_val[i][1] for i in range(df_fits.shape[1])]).T
        # get weighted objective function values based on max/min
        df_fits = df_fits * (weights.values ** [i[0] for i in self.fit_val])
        # scale all values in range 1-2
        df2 = [((df_fits[i] - df_fits[i].min()) / (df_fits[i].max() - df_fits[i].min())) + 1 for i in range(df_fits.shape[1])]
        # inverse min columns
        df2 = pd.DataFrame([df2[i]**self.fit_val[i][0] for i in range(len(df2))]).T
        # rescale all values in range 1-2
        df2 = pd.DataFrame([((df2[i] - df2[i].min()) / (df2[i].max() - df2[i].min())) + 1 for i in range(df2.shape[1])])
        
        fitnesses = list(df2.sum())

        if choice == "Roulette":
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
        else:
            fits_sort = sorted(fitnesses, reverse=True)
            best = [deepcopy(population[fitnesses.index(fits_sort[i])]) for i in range(min(num, len(population)))]
            return best
    
    def search(self, n_generations=20, init_ratio = 0.35, crossover_ratio = 0.35):
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
        def fit_eval(invalid_ind, fit_list):
            epop, fit_list = [i[0] for i in fit_list], list(fit_list)
            new_pop = []
            if invalid_ind: 
                invalid_ind = [i for i in invalid_ind if i not in epop]
                obval = [self.pre_eval(i) for i in invalid_ind]
                rev_smiles, fitnesses = [i[0] for i in obval], [i[1] for i in obval]
                for ind, fit, r_smi in zip(invalid_ind, fitnesses, rev_smiles):
                    if fit is not None:
                        fit_list.append((deepcopy(ind), fit, r_smi))
                        new_pop.append(deepcopy(ind))
            return tuple(fit_list), new_pop

        if init_ratio >=1 or crossover_ratio >=1 or (init_ratio+crossover_ratio)>=1: raise Exception("Sum of parameters init_ratio and crossover_ratio should be in the range (0,1)")
        if self.population is not None:
            pop = self.population
            fit_list = self.fit_list
        else:
            pop = self.pop_generator(n=self.pop_size)       # list of tuples
            fit_list = ()
        
        # Evaluate the initial population
        fit_list, pop = fit_eval(pop, fit_list)

        total_pop = []
        for xg in range(n_generations):
            cross_pop, mutant_pop, co_pop, psum = [], [], [], len(fit_list)
            # Generate crossover population
            co_pop = self.select(pop, fit_list, self.crossover_size)
            shflist = pop + total_pop
            random.shuffle(shflist)
            c_total = co_pop + shflist
            for child1, child2 in zip(c_total[::2], c_total[1::2]):
                if (len(fit_list) - psum) >= self.crossover_size: break
                c1, c2 = self.crossover(child1, child2)
                epop = [i[0] for i in fit_list]
                if c1 in epop or c2 in epop or c1==c2: continue
                fit_list, new_cpop = fit_eval([c1, c2], fit_list)
                cross_pop.extend(new_cpop)
            
            # Generate mutation population
            if self.algo == 4:
                mu_pop = self.select(cross_pop, fit_list, self.mutation_size)
            else:
                mu_pop = self.select(pop, fit_list, self.mutation_size)
            
            pre_mu = len(fit_list)
            for mutant in mu_pop + cross_pop + total_pop + pop:
                if (len(fit_list) - pre_mu) >= self.mutation_size: break
                mt = self.custom_mutate(mutant)
                if mt in [i[0] for i in fit_list]: continue
                fit_list, new_mpop = fit_eval([mt], fit_list)
                mutant_pop.extend(new_mpop)
            
            # Select the next generation individuals
            total_pop = pop + cross_pop + mutant_pop
            if self.algo == 2:
                pop = self.select(total_pop, fit_list, self.pop_size)
            elif self.algo == 3:
                p1 = self.select(pop, fit_list, int(init_ratio*self.pop_size), choice="best")
                p2 = self.select(cross_pop, fit_list, int(crossover_ratio*self.pop_size), choice="best")
                p3 = self.select(mutant_pop, fit_list, self.pop_size-len(p1)-len(p2), choice="best")
                pop = p1 + p2 + p3
            else: pop = self.select(total_pop, fit_list, self.pop_size, choice="best")
            
            # Saving each generation
            if xg==0:
                df_gen = pd.DataFrame([i for i in fit_list])
            else:
                df_gen = pd.DataFrame([ i for i in fit_list[-(self.crossover_size+self.mutation_size):]])
            df_gen = df_gen[[2, 1]]
            df_gen.columns = ['Canonical SMILES', 'Fitness Values']
            fname = '/generation_' + str(xg+1) + '.csv'
            df_gen.to_csv(os.path.join(self.output_dir + fname), index=None)

        self.population = pop    # stores best individuals of last generation
        self.fit_list = fit_list

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

def print_ll(sentence, output_dir):
    """Print to logfile.

    Parameters
    ----------
    sentence: str
        string to be printed to logfile

    Returns
    -------

    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    logfile = open(os.path.join(output_dir+'/logfile.txt'),'a')
    print(sentence)
    logfile.write(str(sentence)+"\n")

def print_lle(sentence, output_dir, msg="Aborting the run"):
    """Print to both error file and logfile and then exit code.

    Parameters
    ----------
    sentence: str
        string to be printed to error file and logfile

    Returns
    -------

    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    logfile = open(os.path.join(output_dir+'/logfile.txt'),'a')
    error_file = open(os.path.join(output_dir+'/error_file.txt'),'a')
    print(sentence)
    logfile.write(sentence+"\n")
    error_file.write(sentence+"\n")
    sys.exit(msg)
    