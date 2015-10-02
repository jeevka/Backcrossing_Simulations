
####################################################################################################
# Importing Python Modules
####################################################################################################
from __future__ import division
import sys
import random
import operator
from copy import copy
from copy import deepcopy
import math
import re
import numpy as np
from heapq import heapify
from numpy import *
#####################################################################################################
# Importing user defined modules.
#####################################################################################################
#import Yeast_Simulator
#import Mutations

################################## Intialization ####################################################
# Intialization of the individuals before the cell division starts.
# 1. Initial population contains cell at different ages. So the cells divide at different speeds.
# 2. 50% cells are virgin cells. 25% are one time divided. 12.5% are two times divided. and etc ...
#####################################################################################################
def initialize_population(n_individuals,n_generations,cell_division_max):
    size = n_individuals * 2**n_generations
    genotypes = ones((size,6),dtype=int32)
    
    n_generations = 16

    for n in xrange(n_individuals):
        # Initializing the cell division time and number of generations.
        genotypes[n] = [cell_division_max,cell_division_max,n_generations,0,0,1]
        
    # Making the Point mutation and gene duplication Dict of list 
    PM = {}; 
    for i in range(size):
        PM[i] = []
        
    return genotypes,PM

#############################################################################################
# Group the cells based on their cell division time.
#############################################################################################
def group_cells(genotypes,n_individuals):
    cell_group = {}
    n_individuals = int(n_individuals)
    for i in xrange(n_individuals):
        if cell_group.has_key(genotypes[i][0]):
            cell_group[genotypes[i][0]].append(i) 
        else:
            cell_group[genotypes[i][0]] = [i]
            
    return cell_group

#################################################################################################
# Updating the next cell divison time and age
#################################################################################################
def update_CDT_age(parent,daughter):
    # Updating the next cell division time
    parent[0] = int(parent[0]) + int(parent[1])
    daughter[0] = int(daughter[0]) + int(daughter[1])
    # Updating the age
    parent[2] = parent[2] - 1 
    daughter[2] = daughter[2] - 1

    return parent,daughter

####################################################################################################
# Convert the fitness values to beneficial and deleterious bases on their proportions
####################################################################################################
def make_fitness_proportions(fitness,MP):
    """
    Input:
        List cell division times
        Proportion of beneficial mutations
    
    Output:
        List of fitness values converted to beneficial or deleterious
    """
    import random as r
    n1 = int(len(fitness) * (100-MP)/100)
    n2 = len(fitness)
    pop = xrange(n2)
    sampled = r.sample(pop,n1)
    for i in sampled:
        fitness[i] *= -1

    return fitness

####################################################################################################
# Convert the fitness values to minutes from random gamma values
####################################################################################################
def covert_to_minutes(fitness):
    """
    Input:
        List of random values from gamma distribution
    
    Output:
        gamma distribution values are converted to mins
    """
    fitness_mins = []
    for i in fitness:
        fitness_mins.append(calculate_cell_division_time_1(i))
    
    return fitness_mins

####################################################################################################
# Convert the fitness values to minutes
####################################################################################################
def calculate_cell_division_time_1(fitness):
    """
    Input:
        one random value from Gamma distribution
    
    Output:
        converted value back
    """
    cell_division_max = 180
    cell_division_min = 90
    #change_CDT = (cell_division_max-cell_division_min) * alpha
    
    change_CDT = (fitness/(1+fitness)) *  cell_division_max

    if change_CDT > 90:
        change_CDT = 89
        
    return change_CDT

####################################################################################################
# Values bigger than the difference between WT and Mutatant cant be included due to the design of our model
####################################################################################################
def truncate_fitness_effects(fitness,trunc):
    """
    Input:
        List of fitness values
    
    Output:
        Truncated values
        e.g. all the values >89 will be removed
    """
    new_fitness = []
    for i in fitness:
        if i <= trunc:
            new_fitness.append(i)
    
    return new_fitness

####################################################################################################
# Calculating new cell division times based on current cell division time and fitness
####################################################################################################
def calculate_cell_division_time(current_CDT,alpha):
    """
    Input:
        current cell division time of a cell
        Fitness effect of the new mutation
    
    Output:
        New cell division time calculated based on exisiting cell division time
        and fitness effect of the new mutation
    """
    # Cell Division time to "Alpha"-refer the coding documents.
    cell_division_max = 180
    cell_division_min = 90
    CDT_change = 0
    if alpha > 0:
        CDT_change = (current_CDT - cell_division_min)/(cell_division_max - cell_division_min)
        CDT_change *= alpha
        new_CDT = current_CDT - CDT_change
    else:
        new_CDT = current_CDT + abs(alpha)

    return  new_CDT