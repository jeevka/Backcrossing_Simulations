from __future__ import division
import sys
import random
from copy import deepcopy
from heapq import *
import heapq
from stdlib cimport *
from numpy import *
import numpy as np
cimport numpy as np
cimport cython
import gc
import re  
import Yeast_Simulator
import Sub_Functions
import Mutations
import BackCrossing

cdef int cell_max_age = Yeast_Simulator.cell_max_age
cdef double mutation_rate = Yeast_Simulator.mutation_rate
cdef double beneficial_mutation_rate = Yeast_Simulator.beneficial_mutation_rate
cdef int n_generations = Yeast_Simulator.number_of_generations

# Cython Imports
cdef extern from "stdlib.h":
    long c_libc_random "random"()
    void c_libc_srandom "srandom"(unsigned int seed)

##############################################################################################
# Each cell divides at differnt speed. 
##############################################################################################
@cython.boundscheck(False)
def asymmetrical_cell_division(genotypes,cell_groups,PM_fitness,
                                PM,PM_strand,n_generations,n_BN,n_individuals,n_mut,point_mutation,ploidy,text,div_type,s_size,fix_cal,PM_Track):
    
    # Cython variable Declarations
    cdef int dps,i,time,tp
    cdef int small
    cdef int k,h,n_cells, P_cdt, C_cdt

    import heapq
    groups = cell_groups.keys()
    groups.sort()
    time = heapq.heappop(groups) 
    
    # desired population size.
    dps = n_individuals * 2**n_generations

    # PM_TRACK_FREQ
    PTF = {}
    
    for i in xrange(n_BN):
        
        if i == 0:
            # For Fixation analysis
            PM_File = open("PM.txt","w")
        else:
            PM_File = open("PM.txt","a")
        
        #PM_File = open("PM.txt","w")
        
        #print "Bottleneck :", i
        tp = n_individuals
        # fastest cell division time
        small = genotypes[0][1]
        # slowest cell division time
        big = genotypes[0][1]
        while tp < dps: 
            # Iterating through the cell groups.
            n_cells = len(cell_groups[time])
            for h in xrange(n_cells):
                k = cell_groups[time][h]
                # Checking whether the cell can divide or not based on how many times it has divided and Haploids and diploids.
                divide = decide_on_division(genotypes[k],div_type)
                if divide == 1:
                    # Copying the parents attributes to child.
                    genotypes[tp] = genotypes[k]
                    genotypes[tp][2] = cell_max_age
                    genotypes[tp][5] = genotypes[k][5]
                    
                    PM[tp] = PM[k][:]
                    
                    # Introduce mutations
                    if random.random() <= mutation_rate:
                       # Sending mother and daughter cells to introduce mutations. 
                       genotypes[k],genotypes[tp],PM, PM_fitness, PM_strand = Mutations.introduce_mutation(genotypes[k],genotypes[tp],
                                                              beneficial_mutation_rate,k,tp,point_mutation,n_mut,PM,PM_fitness,PM_strand,ploidy,fix_cal)
                    
                       n_mut += 1
                       
                    # Reduce the Age for the parent cell
                    genotypes[k][2] = reduce_age(genotypes[k][2])
                    
                    # Updating  the next cell division time based on the cell division time
                    # for both mother and daughter cells.
                    (genotypes[k][0],genotypes[tp][0]) = update_next_CDT(genotypes[k][0],genotypes[k][1],genotypes[tp][0],genotypes[tp][1])                  
                    
                    # Group the cells from current time point
                    cell_groups, groups = group_the_cells(cell_groups,groups,genotypes[k][0],genotypes[tp][0],k,tp)
                    
                    tp = cell_count(tp)
                    
                    # Come out of the loop once its reached the desired population size.
                    if tp == dps:
                       break
                #else:
                #    genotypes[k][0] += genotypes[k][1]

            del cell_groups[time]
            time = heapq.heappop(groups)

        # Sampling - Bottleneck
        if i != n_BN-1: 
            # Random smapling            
            (genotypes,cell_groups,PM,PTF) = sampling_individuals(genotypes,n_individuals,i,PM,PM_File,n_BN,PTF)
            
            # calculate_mean_cell_division_time(genotypes,i,n_individuals)
            #PM_File.close()
            #freq_above_0p5,freq_above_0p3,freq_above_0p1 = fixation(i,n_individuals,n_mut,PM_fitness,text)
            
        else:
            (genotypes1,cell_groups,PM,PTF) = sampling_individuals(genotypes,s_size,i,PM,PM_File,n_BN,PTF)
        
            if fix_cal == 1:
                calculate_mean_number_of_genetic_variations(genotypes1,s_size)
                print "No. of Point mutations:",n_mut,len(PM_fitness)
                #calculate_haplotype_freq(genotypes1,PM,n_individuals)
            
            #print "No. of mitoses:",tp
            PM_File.close()
            
            # Fixation analysis: The conditon is to avoid unnecessary fixation calculation where its not needed.
            # Because, Fixation calculation is one of the heavy and time consuming process.
                        
            if fix_cal == 1:
                freq_above_0p5,freq_above_0p3,freq_above_0p1 = fixation(i,n_individuals,n_mut,PM_fitness,text,PTF)
                calculate_haplotype_freq(genotypes1,PM,n_individuals,freq_above_0p5)
            
            return genotypes1,PM,PM_fitness,PM_strand,n_mut

        # Its a small mess up. Have to find a better way.
        groups = []
        groups = cell_groups.keys()
        time = heapq.heappop(groups)

#################################################################################################
# Haplotype Frequency and the CDT for that Haplotype
#################################################################################################
cdef calculate_haplotype_freq(genotypes,PM,n_indi,freq_above_0p5):
    hap_freq = {}; CDT = {};
    
    # Same Haplotype: different CDT
    SH_D_CDT = 0
    
    for i in xrange(n_indi):
        key_text = "" 
        for j in xrange(len(PM[i])):
            if j != len(PM[i])-1:
                key_text += str(PM[i][j]) + "-"
            else:
                key_text += str(PM[i][j])
        CDT_temp = genotypes[i][1]
        if hap_freq.has_key(key_text):
            hap_freq[key_text] += 1
            if CDT_temp != CDT[key_text]:
                SH_D_CDT += 1
        else:
            hap_freq[key_text] = 1
            CDT[key_text] = CDT_temp
            
    print "Number of different Haplotypes:",len(hap_freq)
    print "Same Haplotype with different CDT:",SH_D_CDT
    
    # Number of distinct Haplotypes based on Mutations above 5% after ASE
    # Distinct Haplotypes
    DH = {}
    for i in hap_freq:
        for j in freq_above_0p5:
            if re.search(str(j),i):
                DH[i] = i
    
    print "Number of mutations above frequency 0p5:",len(freq_above_0p5)                
    print "Number of different haplotypes based on above 0p5:",len(DH)            
        
    # Priting the Haplotypes which are having high freq
    for i in hap_freq:
        freq = hap_freq[i]/n_indi
        print "ALL_Haplotype_ASE\t",i,"\t",hap_freq[i]/n_indi,"\t",CDT[i]
            
    return 0
        

#################################################################################################
# Haplotype Frequency and the CDT for that Haplotype
#################################################################################################
def calculate_haplotype_freq_1(genotypes,PM,n_indi,N_BC):
    N_BC += 1
    hap_freq = {}; CDT = {};
    # Same Haplotype: different CDT
    SH_D_CDT = 0
    for i in xrange(n_indi):
        key_text = "" 
        for j in xrange(len(PM[i])):
            if j != len(PM[i])-1:
                key_text += str(PM[i][j]) + "-"
            else:
                key_text += str(PM[i][j])
        CDT_temp = genotypes[i][1]
        if hap_freq.has_key(key_text):
            hap_freq[key_text] += 1
            if CDT_temp != CDT[key_text]:
                SH_D_CDT += 1
        else:
            hap_freq[key_text] = 1
            CDT[key_text] = CDT_temp
            
    #print "Number of different Haplotypes:",len(hap_freq)
    #print "Same Haplotype with different CDT:",SH_D_CDT
        
    k = 0
    for i in hap_freq:
        k += 1
        print "Haplotype\t",k,"\t",hap_freq[i]/n_indi,"\t",CDT[i]
        
    # Priting the Haplotypes which are having high freq
    text = "S_Haplotype_ABC" + str(N_BC)
    text1 = "ALL_Haplotype_ABC" + str(N_BC)
    for i in hap_freq:
        freq = hap_freq[i]/n_indi
        
        print text1,"\t",k,"\t",hap_freq[i]/n_indi,"\t",CDT[i],"\t",i
        
        if freq >= 0.01:
            print text,"\t",k,"\t",hap_freq[i]/n_indi,"\t",CDT[i],"\t",i
            
    return 0

#################################################################################################
# Random Sampling the desired number of individuals after certain number of generations.
#################################################################################################
cdef sampling_individuals(genotypes,int n_individuals,int i,PM1,PM_File,n_BN,PTF):
    # Making the Point mutation and gene duplication Dict of list 
    size = len(genotypes)
    PM = {}; 
    for i1 in range(size):
        PM[i1] = []
    
    # Sampling the number of IDs from genotypes dict.
    cdef int n_cells,j
    cdef float total_time = 0
    size = len(genotypes)
    #sampled_ids = random.sample(genotypes,n_individuals)
    sampled_ids = random.random_integers(size-1,size=n_individuals)
    sampled_group = ones((size,6),dtype=int32)
    cell_mutations_temp = []
    cell_gene_duplications_temp = []
    cell_groups = {}

    # Choosing the sampled individuals.
    n_cells = len(sampled_ids)
    n_alive_cells = 0
    for j in xrange(n_cells):
        k = sampled_ids[j]
        if genotypes[k][2] != 0:
            total_time += genotypes[k][1]
            n_alive_cells += 1
        sampled_group[j] = genotypes[k]
        
        # Copying PM and GD info of cells
        PM[j] = PM1[k][:]
        
        
        # Write it in a file
        #if i == n_BN -1:
        PTF = print_in_file(PM_File,PM[j],i,PTF)

        # Storing the mutations in sampled cells.
        # cell_mutations_temp.append(cell_mutations[k])
        # Storing the gene duplications in sampled cells.
        # cell_gene_duplications_temp.append(cell_gene_duplications[k])
        
        # Group the cells based on their cell division time.
        kl = int(genotypes[k][0])
        if cell_groups.has_key(kl):
            cell_groups[kl].append(j)
        else:
            cell_groups[kl] = [j]
    
    # This is to indicate the end of each bottleneck
    #PM_File.write(str("##\n"))
    
    Mean_CDT = total_time/n_alive_cells
    
    LSC = calculate_LSC(Mean_CDT)
    print "LSC\t",i,"\t",Mean_CDT,"\t",LSC
    
    return sampled_group,cell_groups,PM,PTF

##############################################################################################
# Calculating LSC from Mean Cell divisiobn time after each bottleneck/serial transfer
# In otherwords, LSC is relative growth rate of a mutatnt cell growing in stress condition
# relative to WT cells growing in Normal conditions 
##############################################################################################
cdef calculate_LSC(mean_CDT):
        """
        Input:
            Mean cell division time of a population
        
        Output:
            Growth rate relative to WT cell growing in normal growth media
        """
        
        # For Below formula refer to Jonas Warringer's paper.
        LSC = (math.log(180/60) - math.log(mean_CDT/60))/2
        LSC = LSC/0.3465735902799727        
        
        return LSC
        
##############################################################################################
#
##############################################################################################
cdef print_in_file(FH,PM,BN,PTF):
    if len(PM) > 0 :
        FH.write(str(BN))
        FH.write(",")
        for i in PM:
            FH.write(str(i))
            FH.write(",")
    
            if PTF.has_key(i):
                PTF[i][BN] += 1
            else:
                PTF[i] = [0 for j in xrange(50)]
                PTF[i][BN] = 1
        FH.write("\n")
        
    return PTF

###############################################################################################
# To find mean number of mutations, gene deletions, gene duplications.
###############################################################################################
cdef calculate_mean_number_of_genetic_variations(genotypes,int tp):
    cdef int i,mean_mutation = 0
    cdef int mean_gene_duplication = 0

    for i in xrange(tp):
        mean_mutation += genotypes[i][3]
        mean_gene_duplication += genotypes[i][4]
    
    print "Mean Mutation:",mean_mutation,tp,(mean_mutation)/float(tp)
    
    return

##############################################################################################
# To decide whether the cell can divide based on AGE and DIVISION TYPE(either Hap or Dip)
# Its mainly to manage the EXCLUSIVE DIPLOID and EXCLUSIVE HAPLOID division during Backcrossing
##############################################################################################
cdef decide_on_division(genotypes,div_type):
    if genotypes[2] > 0:
        if div_type == 1:
            return 1
        else:
           if genotypes[5] == 2:
              return 1
           else:
              return 0
    else:
        return 0
            
def TRACK_ALL_MUTATIONS(mut_0p01,PTF,n):
    for i in mut_0p01:
        if PTF.has_key(i):
            for j in xrange(50):
                if PTF[i][j] != 0:
                    print "Fixation\t",i,"\t",j,"\t",PTF[i][j]/n
                    
    return 0
    
###############################################################################################
# To calculate the Fixaed number of mutations and gene duplications.
# Here we mainly split,extract and format the data
# By calling the cal_fix we calculate the rate of fixation
###############################################################################################
cdef fixation(BN,n,n_mut,PM_fitness,text,PTF):
    BN += 1
    
    # Read the Mutation output File
    PM = read_data() 
    
    # Grep Mutations from last Bottleneck
    PM_fix = greping(BN,PM)

    # Split the mutations and store them in list
    Last_BN_PM = splitting(PM_fix)
    
    # Take the unique mutation numbers
    U_PM = unique_grep(Last_BN_PM)
    
    # Separate the mutations 
    freq_above_0p5,freq_above_0p3,freq_above_0p1 = cal_fix(Last_BN_PM,U_PM,n,n_mut,"PM",PM_fitness,text)

    #time_to_fix(PM,BN,n_mut)
    
    #TRACK_ALL_MUTATIONS(freq_above_0p1,PTF,n)
    
    return freq_above_0p5,freq_above_0p3,freq_above_0p1


cdef time_to_fix(data,BN,n_mut):
     # Part I: Seperate mutations based on Bottlenecks.
     BN_Mut = []
     U_BN_Mut = []
     for i in xrange(1,BN+1):
        # Grep the data for particular Bottleneck
        temp1 = greping(i,data)
        
        # Splitting the data
        temp2 = splitting(temp1)
        
        # Whole mutation set
        BN_Mut.append(temp2)

        # Unique mutations
        temp3 = unique_grep(temp2)
        U_BN_Mut.append(temp3)
      
     
     """   
     # Part II: Check Each mutation for introduction and extinction time
     for i in range(n_mut):
         t2 = 0; t1 = 0
         for j in range(BN):
             t2 = 0
             for k in range(len(U_BN_Mut[j])):
                 if U_BN_Mut[j][k] == i and t1 == 0:
                    start = j
                    t1 = 1
                 if U_BN_Mut[j][k] == i:
                    t2 = 1
                    
             if t1 == 1 and t2 == 0:
                end = j; time_extinct = (end - start)*5
                print i,"\t",start,"\t",end,"\t",time_extinct
                t2 = 1
                break 
             
         if t1 == 1 and t2 == 1 and j == BN-1:
            end = j;time_extinct = (end - start)*5
            print i,"\t",start,"\t",end,"\t",time_extinct
         
         if t1 == 0 and t2 == 0:
            print  i,"\t0\t0\t0"
     """   
               
###############################################################################################
# To read mutation and gene duplication data
###############################################################################################
cdef read_data():    
    PM1 = open("PM.txt","r")
    PM = PM1.readlines()
    PM1.close()

    return PM
    
cdef Track_The_Mutations(mut):
    data = read_data()
    MUT_FREQ= {i:[0] for i in xrange(50)}
    for i in mut:
       for j in data:
        temp = j.split(",")
        BN = int(temp[0])
        temp2 = temp[1:len(temp)-1]
        if i in temp2:
            MUT_FREQ[i][BN] += 1
        
    print MUT_FREQ
    sys.exit()
    

###############################################################################################
# To calculate the Fixated number of mutations.
###############################################################################################
cdef cal_fix(PM,U_PM,n,n_mut,string,fitness,text):
    freq_above_0p5 = {}
    freq_above_0p3 = {}
    freq_above_0p1 = {}
    n_fix = 0 
    for i in U_PM:
         m = 0 
         for j in PM:
            if i == j:
               m += 1
         fix = float(m)/float(n)
        
         if fix >=0.01:
            freq_above_0p1[i] = fix
            
         if fix >=0.03:
            freq_above_0p3[i] = fix
            
         if fix >=0.05:
            freq_above_0p5[i] = fix
            
         print "PM\t",i,"\t",fix,"\t",fitness[i]         
         m = 0
         
         # Check how many muations are fixed
         if fix == 1:
            n_fix += 1
    txt = "No. of fixated " + str(string) + ":"        
    print txt,n_fix
    if n_fix != 0:
        print "Rate of fixation:",float(n_fix)/float(n_mut)
    
    
    return freq_above_0p5,freq_above_0p3,freq_above_0p1   

###############################################################################################
# To extract the unique mutations from the last bottleneck.
###############################################################################################    
cdef unique_grep(PM):
    set = {} 
    map(set.__setitem__, PM, []) 
    
    return set.keys()

###############################################################################################
# To split the mutation numbers and store them in a list from the last bottleneck
###############################################################################################    
cdef splitting(PM):
    Last_BN_PM = []
    for j in PM:       
        temp1 = re.split(",",j)
        for k in xrange(1,len(temp1)-1):    
            Last_BN_PM.append(int(temp1[k]))

    return Last_BN_PM
    
###############################################################################################
# To grep the mutation profiles from the last bottleneck
###############################################################################################    
cdef greping(BN,PM):
    BN -= 1
    Last_BN = []
    for i in xrange(len(PM)):
        pattern = "^" + str(BN) + ","
        if re.search(pattern,PM[i]):
           Last_BN.append(PM[i].strip())
    
    return Last_BN
           
###############################################################################################
# To update the next cell division time.
###############################################################################################
cdef inline update_next_CDT(int a1, int a2, int b1, int b2):
     a1 = a1 + a2
     b1 = b1 + b2
     return a1,b1
     
cdef inline reduce_age(int a):
     a = a - 1
     return a
     
cdef inline cell_count(int tp):
     tp = tp + 1
     return tp

###############################################################################################
# To group the cels
###############################################################################################
def group_the_cells(cell_groups,groups,P_cdt,C_cdt,k,tp):
    if P_cdt == C_cdt:               
        # Grouping the mother cells
        if cell_groups.has_key(P_cdt):
          cell_groups[P_cdt].append(k)
          cell_groups[C_cdt].append(tp)
        else:
          cell_groups[P_cdt] = [k]
          cell_groups[C_cdt].append(tp)
          heapq.heappush(groups,P_cdt)
    else:
        # Grouping the mother cells
        if cell_groups.has_key(P_cdt):
          cell_groups[P_cdt].append(k)
        else:
          cell_groups[P_cdt] = [k]
          heapq.heappush(groups,P_cdt)                    
        # Grouping the daughter cells.
        if cell_groups.has_key(C_cdt):
          cell_groups[C_cdt].append(tp)
        else:
          cell_groups[C_cdt] = [tp]
          heapq.heappush(groups,C_cdt)
    
    return cell_groups,groups