import pandas as pd
import numpy as np
from collections import defaultdict


from .multiply import bool_multiply

class Multi_value_minterm_on:
    
    def __init__(self, minterm, coverage):
        self.minterm = tuple(x for x in minterm)
        self.coverage = frozenset(x for x in coverage)
        
    def __str__(self):
        return('{0}, {1}'.format(
                                    self.minterm, 
                                    self.coverage
                                    
                                    ))
    def __repr__(self):
        return str(self)
        
    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (self.minterm == other.minterm and 
                self.coverage == other.coverage) 
                
    def __hash__(self):
        return hash((self.minterm, self.coverage))
        
    def can_be_reduced(self,other):
        n = len(self.minterm)
        diff = 0
        for i in range(0,n):
            if (self.minterm[i] != other.minterm[i] and 
                 (self.minterm[i]== int(-1) or other.minterm[i] == int(-1))):
                diff = diff+1
        if diff == 1:
            return(True)
        else:
            return(False)
      
    def reduce(self, other):
        if not self.can_be_reduced(other):
            return None
  
        n = len(self.minterm)
        new_minterm = list()
        for i in range(0,n):
            if(self.minterm[i] == other.minterm[i]):
                new_minterm.append(self.minterm[i])
            else:
                new_minterm.append( int(-1))
        return Multi_value_minterm_on(new_minterm, 
                       self.coverage.union(other.coverage))

        
    
    def reduce_with_off_set(self, off_matrix):
        res = []       
        for ind,elm in enumerate(off_matrix.data):

        
            new_minterm = list()
            n= len(self.minterm)

            for i in range(0,n):
               
                if(self.minterm[i] != elm[i]):
                    new_minterm.append( self.minterm[i])
                else:
                    new_minterm.append(int(-1))
            res.append(Multi_value_minterm_on(new_minterm, 
                       self.coverage))

        return res


class On_off_reduction_matrix:
    
    def __init__(self,matrix):
        self.matrix = matrix
    
    
    def __str__(self):
        return('{0}'.format(self.matrix))
                                    
    def __repr__(self):
        return str(self)
    
    def reduction(self):
        any_reduction = True
      
        while(any_reduction):
            any_reduction = False
            used = [False] * len(self.matrix)
            new_matrix = []
            for i in range(len(self.matrix)):
                for j in range(i+1, len(self.matrix)):
                    if used[i] or used[j]:
                        continue
                    m1 = self.matrix[i]
                    m2 = self.matrix[j]
                    if m1.can_be_reduced(m2):
                        new_matrix.append(m1.reduce(m2))
                        any_reduction = True
                        used[i] = True
                        used[j] = True
                        
            for is_used, m in zip(used, self.matrix):
                if not is_used:
                    new_matrix.append(m)
            self.matrix = new_matrix
                    
        return self.matrix
            
        
         
        

class Multi_value_off_matrix:
    
    def __init__(self, data):
        self.data = data
        
       
    def __len__(self):
        return len(self.data)
    def __str__(self):
        return('{0}'.format(self.data))
                                    
                                    
    def __repr__(self):
        return str(self)
        
   
        
   
def set_labels(term_object,labels):
    implicant = []
    for i,j in zip(term_object[0],labels):
        if i == int(-1):
            continue
        else:
            if i == 1:
                implicant.append(str(j.upper())) 
            elif i == 0:
                implicant.append(str(j.lower())) 
    implicant_complete = '*'.join(implicant)
    return (implicant_complete, term_object[1])



def on_off_grouping(table, output, multi_output = False):

    data_positive =table[table[output] == 1]
    inputs =  [x for x in table.columns if x!= output]
    data_on =  np.array(data_positive[inputs])
    ind_on = [int(x) for x in table[table[output] == 1].index]
    data_negative = table[table[output] == 0]
    
    data_off = np.array(data_negative[inputs])
    
    onset = [Multi_value_minterm_on(row,{ind}) for row,ind in
             zip(data_on,ind_on)]
    offset = Multi_value_off_matrix(data_off)
      
    return onset, offset

def prepare_m_for_bool_multiply(m_in):
    if len(m_in) == 0:
        raise RuntimeError()
    coverage = m_in[0].coverage
    m = [x.minterm for x in m_in]
    return coverage, m
    

def reduction(onset, offset):
    on_off_matrixes = []
    for minterm in onset:
        on_off_matrixes.append(On_off_reduction_matrix(
    
                                    minterm.reduce_with_off_set(offset)))        
    
    imp_cov_dict = defaultdict(lambda: set())
    
    for matrix in on_off_matrixes:
        reduced_matrix = matrix.reduction()
 
        coverage,m = prepare_m_for_bool_multiply(reduced_matrix)
        res = bool_multiply(m)
        
        
        for imp in res:
            imp_cov_dict[imp].update(coverage)
    return imp_cov_dict
   