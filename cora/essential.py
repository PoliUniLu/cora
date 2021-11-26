#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from .multiply import bool_multiply
from functools import reduce



def create_onset_groups(onset):
    max_len = len(onset[0])+1
    groups = list()
    for i in range(0,max_len):
        groups.append([])
    for x in onset:
        groups[sum(1 for elm in x if elm!=0)].append(x)
    return groups

def create_offset_groups(offset):
    max_len = len(offset[0])+2
    groups = list()
    for i in range(0,max_len):
        groups.append([])
    for x in offset:
        groups[sum(1 for elm in x if elm!=0)].append(x)
    return groups

def are_reducable(el1,el2):
    sum = 0
    for i,j in zip(el1,el2):
            if i!=j:
                sum +=1
    return sum==1
def reduction(el1,el2):
    e = [-1]*len(el1)
    for ind,(x,y) in enumerate(zip(el1,el2)):
            if x!=y:
                e[ind] = x
    return e
            
# compare each element from P(x) with S(x-1) and S(x+1)
def compare_the_groups(P,S,multi_value = False):
     l = len(P)
     res = []
     for i in range(0,l):
         res.append([])
         
     for ind in range(0,l):
         for ind2,el1 in enumerate(P[ind]):
             res[ind].append(set())
             if S[ind+1]:
                 for el2 in S[ind+1]:
                    
                      if are_reducable(el1,el2):
                       
                              
                              res[ind][ind2].add(tuple(reduction(el1,el2)))
                              
             if multi_value:
                if S[ind]:
                    for el2 in S[ind]:
                        if are_reducable(el1,el2):
                              res[ind][ind2].add(tuple(reduction(el1,el2)))
                       
                              
             if S[ind-1]:
                 for el2 in S[ind-1]:
                     if are_reducable(el1,el2):
                              res[ind][ind2].add(tuple(reduction(el1,el2)))
     return res
                          
def multiply_out(groups):
    res = []
    for group in groups:
        if group:
            for impl_subgroup in group:
               res.append( bool_multiply(list(impl_subgroup)))
                
    return res
                          


def on_off_grouping_mo(data,output):
   
    onset_tmp = data[data[output] == 1]
    onset = np.array(onset_tmp[[x for x in data.columns if x!=output]])
    
    offset_tmp = data[data[output]==0]
    offset = np.array(offset_tmp[[x for x in data.columns if x != output]])
    
    return onset, offset    

def merge_elements(el1, el2):
    return tuple(max(x,y) for x,y in zip(el1,el2))

def reduce_group(group):
    return {reduce(merge_elements, subgroup) for subgroup in group if len(subgroup) > 0}

def combine_an_implicant(groups):
    res = set()
    for group in groups:
        res.update(reduce_group(group))
    return res

def included_in_offset(implicant,data,output):
    offset = data[data[output]==0]
    return any(offset.apply(
               lambda row_series: True if all(y == x for x,y in
                                          zip(row_series.values,
                                              implicant ) if y >=0) 
                                   else False, axis = 1))
               
def reduce_the_onset(essentials,data, output):
    onset = data[data[output] == 1]  
    cov_list = []
    for imp in essentials:
       cov_list.append( onset[onset.apply(
            lambda row_series: True if all(y == x for x,y in 
                                           zip(row_series.values,
                                               imp) if y>=0)
                                            else False,axis = 1 )].index) 
    return cov_list

def transform_to_raw_impl(impl, levels):
    res = [frozenset(range(i)) for i in levels]
    for ind,x in enumerate(impl):
        if x == -1:
            continue
        else:
            res[ind] = frozenset([x])
    return tuple(res)



def get_essential_implicants(data, output,levels):
    onset, offset = on_off_grouping_mo(data,output)
    P = create_onset_groups(onset)
    S = create_offset_groups(offset)
    multi_value = True if max(levels)>2 else False
    groups = compare_the_groups(P,S,multi_value)
    elimination = combine_an_implicant(groups) 
    imps =[imp for imp in elimination if not 
          included_in_offset(imp,data,output)]
    return imps


    
