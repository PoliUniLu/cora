from collections import defaultdict
import pandas as pd
import numpy as np

from .on_off_alg import MultiValueMintermOn, OnOffReductionMatrix, _bool_multiply
from .multiply import _transform_to_raw_implicant

    
class MultiValueMintermOnMo(MultiValueMintermOn):
    
    def __init__(self, minterm, coverage, tag):
        super().__init__(minterm, coverage)
        self.tag = tag
   
    def __str__(self):
        return('{0}, {1},{2}'.format(
                                    self.minterm, 
                                    self.coverage,
                                    self.tag
                                    
                                    ))
    def __repr__(self):
        return str(self)
    
    
    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (self.minterm == other.minterm and 
                self.coverage == other.coverage and
                self.tag == other.tag) 
                
    def __hash__(self):
        return hash((self.minterm, self.coverage, self.tag))
    
     
    def _can_be_reduced(self,other):
        return self.tag == other.tag and super()._can_be_reduced(other)
    
    def _reduce(self, other):
        tmp = super()._reduce(other)
        return MultiValueMintermOnMo(tmp.minterm, tmp.coverage, self.tag)
     
    def reduce_with_off_set(self, off_matrix):
        res = []       
        for ind,(elm,tag) in enumerate(zip(off_matrix[0], off_matrix[1])):
        
            new_minterm = list()
            n= len(self.minterm)

            for i in range(0,n):
               
                if(self.minterm[i] != elm[i]):
                    new_minterm.append( self.minterm[i])
                else:
                    new_minterm.append(int(-1))
                
            new_tag = [x*y for x,y in zip(tag, self.tag)]
           
            if(any(new_tag)):
                
                res.append(MultiValueMintermOnMo(new_minterm, 
                       self.coverage, new_tag))

        return res
    
    
class MultiValueOffMatrixMo:
    
    def __init__(self, data, tags):
        self.data = data
        self.tags = tags
        
       
    def __len__(self):
        return len(self.data)
    def __str__(self):
        return('{0}'.format(self.data))
                                    
                                    
    def __repr__(self):
        return str(self)
    
class OnOffReductionMatrixMo(OnOffReductionMatrix):
    
    def __init__(self, matrix):
        super().__init__(matrix)
        

    
def on_off_grouping_mo(table,outputs):
    data_output = table[outputs]
    onset_index = data_output.apply(lambda row_series:
                                    any(x for x in row_series),axis = 1)
    input_columns = [x for x in table.columns if  x not in outputs]
    inputs = table[input_columns]
    onset_table = inputs[onset_index].apply(lambda row : tuple(int(x) if not np.isnan(x) else None for x in row), axis=1)
    tags_tuples = table[outputs].apply(tuple, axis=1)
    ind_on = [int(x) for x in table[onset_index].index]
    
    
    onset = [MultiValueMintermOnMo(row,{ind},tags) for row,ind,tags in
             zip(onset_table,ind_on,tags_tuples[onset_index])]

    offset_index =  data_output.apply(lambda row_series:
                                     not all(x for x in row_series),axis = 1)
    offset_data = data_output

    offset = table[offset_index]
    offset_data = offset[input_columns].apply(tuple, axis=1)
    offset_tags_tmp = offset[outputs].applymap(lambda x: 1 if x == 0 else 0)
        
    offset_tags = offset_tags_tmp.apply(tuple, axis = 1)

    #offset_tags = offset_tags_tmp.apply()
    offset = offset_data, offset_tags

      
    return onset, offset

def number_to_tuple(x: int, l: int):
    res = []
    for _ in range(l):
        res.append(x & 1 == 1)
        x >>= 1
    return tuple(res)
        

def reduction_mo(onset, offset):
    
    on_off_matrices = []

    for minterm in onset:
        m_res = minterm.reduce_with_off_set(offset)

        m = OnOffReductionMatrixMo(m_res)

        on_off_matrices.append(m)  
   
    

    tag_length = len(next(x for x in onset).tag)
    
    
    tag_imp_dict = defaultdict(list)
    for m in on_off_matrices:
        m_reduced = m.reduction()
        
        for current_tag in map(lambda x: number_to_tuple(x, tag_length), range(1, 2**tag_length)):
            current_minterms = [mt for mt in m_reduced
                                if any(x > 0 and y > 0 for x,y in zip(current_tag, mt.tag))]
            
   
            b_m = _bool_multiply([x.minterm for x in current_minterms])
            tag_imp_dict[tuple([x for x in reversed(current_tag)])].extend(b_m)
            
       
    return tag_imp_dict
   
            
    
        

