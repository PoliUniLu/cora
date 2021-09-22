

import pandas as pd
import numpy as np


class ValuedVariable:
    
    def __init__(self, identificator, value):
        self._ident = identificator
        self._val = value
        
    def __hash__(self):
        return hash((self._ident, self._val))
    
    def __eq__(self, other):
        return self._ident == other._ident and self._val == other._val
    
    def __str__(self):
        return '{}({})'.format(self._ident, self._val)
    
    def __repr__(self):
        return str(self)



def get_levels(data, outputs,inputs):
      inputs = data.drop(outputs,axis=1)
      if len(inputs) == 1:
          dim_corrected = [inputs.iloc[:,0].values]
      else:
          dim = inputs.apply(lambda x: pd.unique(x).tolist(),axis=0,
                             result_type='reduce').array
          
          dim_corrected = []
          for ar in dim:
              if len(ar) > 1:
                  dim_corrected.append(ar)
              else:
                  dim_corrected.append([0,1])
      levels =  [len(x) for x in dim_corrected]
      return levels

def initialze_minterms(levels):
    l= len(levels)
    res = []
    for ind,el in enumerate(levels):
        
        for i in range(el):
            tmp=[-1] * l
            tmp[ind]=i
            res.append(tmp)

    return res
#combination - tuple (-1, 1,-1,-1)
              
def is_an_impicant(data,output,inputs,combination):
     tmp_onset =data[data[output] == 1]
     onset = tmp_onset[inputs]
     tmp_offset = data[data[output] == 0]
     offset = tmp_offset[inputs]

     return (any(onset.apply(
               lambda row_series: True if all(y == x for x,y in
                                          zip(row_series.values,
                                              combination ) if y >=0) 
                                   else False, axis = 1))
               and 
        
            all(offset.apply(
                lambda row_series: False if all(y == x for x,y in
                                          zip(row_series.values,
                                             combination) if y >= 0) 
                                   else True, axis = 1))


            )
from collections import deque


def minterm_value(m1):
    s1 = 0
    for ind,el in enumerate(m1):
        if el != -1:
            s1=s1+ind+1
    return s1       
       
def merge_minterms(m1,m2):
    res = []
    for ind,el in enumerate(m1):
        if el == -1:
            res.append(m2[ind])
        else:
            if m2[ind] == -1:
                res.append(m1[ind])
            else:
                res.append(m1[ind]+m2[ind])
    return res

def extend_minterm(m1,levels):
    res = []
    # Find last defined variable (!= -1)
    ldv = max(i for i,v in enumerate(m1) if v > -1)
    # We need to switch variables starting from ldv+1 
    l = len(m1)
    extentions = []
            
 
    for i  in range(ldv+1, len(m1)):
        for j in range(levels[i]):
            extended_minterm = m1[0:i] + [int(j)] + m1[i+1:l] 
            extentions.append(extended_minterm)
        #m_positive = m1[0:i] + [1] + m1[i+1:l]
        #m_negative = m1[0:i] + [0] + m1[i+1:l]
        #res.extend([m_positive, m_negative])
        res.extend(extentions)
    return res
        
def find_implicants_cubes(data,outputs,inputs):
    levels = get_levels(data,outputs,inputs)
    minterms_inputs = initialze_minterms(levels)
    q = deque(minterms_inputs)
    res = []
    while len(q) >  0:
        x = q.popleft()
        if any(implicant_includes(imp,x) for imp in res):
            
            continue
        elif is_an_impicant(data,outputs,inputs,x):
            res.append(x)       
        else:
            q.extend(extend_minterm(x,levels))

    return res

    
def implicant_includes(x,y):
    # Implicant x includes y
    return all(i == j or i == -1 for i,j in zip(x,y))



def transform_to_raw_imp(impl, levels):
    s = {(ValuedVariable(i, v) 
                        ) for i,v in enumerate(impl)
                        if v != -1}
    res = [frozenset(range(i)) for i in levels]   
    for x in s:
        res[x._ident] = frozenset([x._val])
    return tuple(res)

def minterm_to_str(minterm, levels, labels, tag,multi_output):
    is_multi_level = any(x > 2 for x in levels)
    tmp = [set_to_str(x, y, z, is_multi_level) for x,y,z in
           zip(minterm, levels, labels)]
    res = '{}'.format('*'.join(x for x in tmp if x != ''))
    
    return res if res != '' else '1'

def set_to_str(s,levels,label, is_multi_level):
    if len(s) == levels:
        return ''

    if not is_multi_level:
        if 0 in s:
            return label.lower()
        else:
            return label.upper()
    if len(s) == 1:
        return '{}{{{}}}'.format(label,tuple(s)[0])
    return '{}{{{}}}'.format(label,','.join(str(x) for x in s))

