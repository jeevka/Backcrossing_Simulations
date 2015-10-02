from __future__ import division
import random
import sys
import math
import numpy
import scipy

###################################################################################################
# Import user defined Modules.
###################################################################################################
import Yeast_Simulator
import Sub_Functions as SF
import Cell_Division
import Artificial_Cells
###################################################################################################
# Global variables - Importing from main Modules.
###################################################################################################
# Gamma distribution values for Mutations. Refer the coding document.Or see the refernce Below.

fitness_affecting_mutation_rate = Yeast_Simulator.fitness_affecting_mutation_rate
mutation_fitness = Yeast_Simulator.mutation_fitness_random_fitness

###################################################################################################
# Introducing Mutations Based on Yeast Mutation rate.
###################################################################################################
def introduce_mutation(mother_cell,daughter_cell,beneficial_mutation_rate,k,tp,point_mutation,n_mut,PM,PM_fitness,PM_strand,ploidy,fix_cal):
    beneficial = 0
    deleterious = 0
    fitness_proportion = 0
    
    # Decide whether its fitness altering mutation or neutral.
    # Beneficial Mutation rate is 575 out of 10000 fitness altering mutations.
    # Its 0.088 or 88/1000. Plz refer the coding value documents.
    # From Coding documents: 1.53% of all mutations are Fitness altering mutations.

    if random.randrange(100) < fitness_affecting_mutation_rate:
        # Decide Whether its Adavantageous or Disadvantageous or neutral mutation.
        if random.random() <= beneficial_mutation_rate:
            beneficial = 1
        else:
            deleterious = 1
    
    # fix_cal decides whether GDs carry any fitness effects.
    if fix_cal == 0:
        fitness = 0
    else:
        if random.randrange(100) < fitness_affecting_mutation_rate:        
            fitness = mutation_fitness[n_mut]
        else:
            fitness = 0 
        
    #haplotype = PM_haplotype[n_mut]
    
    # Choose the cell:Parent or child to have mutation.
    if random.random() <= 0.5:
        # Assigning the Chromosome and the position.
        point_mutation.write(str(assign_haplotype(mother_cell[5])))
        point_mutation.write("\t")
        if fitness != 0:
            mother_cell[1] = SF.calculate_cell_division_time(mother_cell[1],fitness)
            
        PM_fitness[n_mut] = fitness
        PM_strand[n_mut] = random.randrange(1,ploidy+1)
        point_mutation.write("\n")
        mother_cell[3] += 1
        
        # Appending the mutation number
        PM[k].append(n_mut)
        
    else:
        point_mutation.write(str(assign_haplotype(daughter_cell[5])))
        point_mutation.write("\t")        
        if fitness != 0:
            daughter_cell[1] = SF.calculate_cell_division_time(daughter_cell[1],fitness)
        PM_fitness[n_mut] = fitness        
        PM_strand[n_mut] = random.randrange(1,ploidy+1)
        
        point_mutation.write("\n")
        daughter_cell[3] +=  1
       
        # Appending the mutation number
        PM[tp].append(n_mut)
    
    beneficial = 0
    deleterious = 0
    
    return mother_cell,daughter_cell,PM,PM_fitness,PM_strand
    
###################################################################################################
# Assigning the Haplotype structure.
###################################################################################################
def assign_haplotype(cell_type):
        
    # Chromosome numbers and the length.
    chrom = {
             1:230208,2:813178,3:316617,4:1531918,
             5:576869,6:270148,7:1090946,8:562643,
             9:439885,10:745745,11:666454,12:1078175,
             13:924429,14:784333,15:1091289,16:948062
             }

    haplotype_structure = {}

    # Choose a chromosome randomly.
    chr_no = random.randint(1,16)
    
    """
    # Check for Haploid or diploid
    if cell_type == 1:
        pair_no = 1
    else:
        if random.random() <= 0.5:
            pair_no = 1
        else:
            pair_no = 2        
    """
    
    # Choose a specific position in chromosome randomly.
    bp_position = random.randint(1,chrom[chr_no])
    haplotype_structure[chr_no] = bp_position
    
    return haplotype_structure