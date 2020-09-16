import numpy as np
import pandas as pd
import itertools
from .petric import find_irredundant_sums,boolean_multiply
from functools import reduce
import string 



COLUMN_LABELS = list(string.ascii_uppercase) + ["AA", "BB", "CC", "DD", "EE", "FF"]

def concatenate_strings(arr):
    return ','.join([str(x) for x in arr])

def inclusion_score(arr):
    return sum(arr)/len(arr)



def count_non_zeros(minterm):
    n=len(minterm)
    non_zeros=0
    for i in range(0,n):
        if minterm[i]=={0}:
            continue
        else:
            non_zeros=non_zeros+1
    return non_zeros

def vector_to_minterm(v):
    return tuple(frozenset({x}) for x in v)

def preprocess_input(table):
    return [vector_to_minterm(x) for x in table]


def create_care_translation_dict(cares, labels):
  return {k:v for k,v in zip(cares, labels)}



def check_element_coverage(row, element):
    for r,e in zip(row, element[0]):
        if r not in e:
            return False
    return True


def is_minterm_subset(m1,m2):
    for x,y in zip(m1[0], m2[0]):
        if not x.issubset(y):
            return False
    return True


def create_groups(table,column_number, cares, outputcolumns, multi_output):
    res=[]
    #column_number define the number of possible nonzero digits
    for i in range(0,column_number+1):
        res.append([])
    index=0
    if not multi_output:
      for row in table:
          non_zeros=count_non_zeros(row)
          res[non_zeros].append(Multi_value_item(row,set([index]) if index in cares else set()))
          index=index+1
      return res
    else:
      nr_of_outputs=len(outputcolumns)
      care_index = 0

      for row in table:
          non_zeros=count_non_zeros(row)
          if index in cares:
              tag={i for i in range(1,nr_of_outputs+1) if outputcolumns[i-1][care_index]==1}
              care_index=care_index+1
          else:
              tag={i for i in range(1,nr_of_outputs+1)}
          if len(tag) == 0:
              index=index+1
              continue
          res[non_zeros].append(Multiple_output_item(row,set([index]) if index in cares else set(),tag))
          index=index+1
      return res




def reduction_step(groups,n,multi_output):
    was_any_reduction=False
    new_groups = []
    for i in range(n+1):
        new_groups.append(set())
    

    for i in range(0,n):
        if len(groups[i])==0:
            continue
        for element1 in groups[i]:
            for element2 in groups[i+1]:
                if element1.can_be_reduced(element2):
                    reduced_element = element1.reduce(element2)
                    new_groups[i].add(reduced_element)
                    was_any_reduction=True
     
    for i in range(0, n+1):
        for element1 in groups[i]:
            for element2 in groups[i]:
                if element1.can_be_reduced(element2):
                    reduced_element = element1.reduce(element2)
                    new_groups[i].add(reduced_element)
                    was_any_reduction=True                
    
    final_implicants=set([])
    for i in range(0,n+1):
        for element1 in groups[i]:
            if(element1.is_reduced is False and (len(element1.coverage) > 0)):
                if multi_output:
                  final_implicants.add((tuple(element1.minterm),frozenset(element1.coverage), frozenset(element1.tag)))
                else:
                  final_implicants.add((tuple(element1.minterm),frozenset(element1.coverage))) 


            
    return({'groups':new_groups,'implicants':final_implicants,'reduction':was_any_reduction})

#--------------------extra_elimination------------------------#
# for complex PI, which were just partially reduced to create several new PIs with totally reduced variables

def extend_with_first(first, arr,multi_output):
  if multi_output:
    return [(tuple([first])+x[0], x[1],x[2]) for x in arr]
  return  [(tuple([first])+x[0], x[1]) for x in arr]



def decomposition(element, levels,multi_output):
    if len(element[0])==0:
        return [element]
    if multi_output:
      tmp = decomposition((element[0][1:], element[1],element[2]), levels[1:],multi_output)
    else:
      tmp = decomposition((element[0][1:], element[1]), levels[1:],multi_output)

    if len(element[0][0]) == 1 or len(element[0][0])==levels[0]:
        return extend_with_first(element[0][0], tmp,multi_output)
    else:
        res = []
        for x in element[0][0]:
            res.extend(extend_with_first({x},tmp,multi_output))
        return res

def fix_coverage_after_decomposition(table, elements,multi_output):
    new_elements = []
    for e in elements:
        new_coverage = frozenset([x for x in e[1] if check_element_coverage(table[x], e)])

        if len(new_coverage) > 0:
          if multi_output:
            new_elements.append(tuple([e[0],new_coverage,e[2]]))
          else:
            new_elements.append(tuple([e[0],new_coverage]))

    return new_elements

#-------------full eliination---------------
#usning the reduction step while there is something to minimize 
#and reducing complex PIs to singular elements

def eliminate_minterms(table,elements, levels,multi_output):
    decomposed = []
    for x in elements:
        decomposed.extend(decomposition(x, levels,multi_output))

    decomposed = fix_coverage_after_decomposition(table, decomposed,multi_output)
    n = len(decomposed)
    was_eliminated=[False]*n
    if multi_output:
      for i in range(n):
          for j in range(i+1,n):
              if is_minterm_subset(decomposed[i],decomposed[j]) and (decomposed[i][2]).issubset(decomposed[j][2]):
                  was_eliminated[i] = True
              elif is_minterm_subset(decomposed[j],decomposed[i]) and (decomposed[j][2]).issubset(decomposed[i][2]):
                  was_eliminated[j] = True
      return [x for we,x in zip(was_eliminated,decomposed) if not we]   
    
    for i in range(n):
        for j in range(i+1,n):

            if is_minterm_subset(decomposed[i],decomposed[j]):
                was_eliminated[i] = True
            elif is_minterm_subset(decomposed[j],decomposed[i]):
                was_eliminated[j] = True
    return [x for we,x in zip(was_eliminated,decomposed) if not we]


    
def find_index2(arr, x):
    return np.where((arr == x).all(axis=1))[0]  



class Chart:

  def __init__(self,data,output_labels,input_labels=None,case_col=None,n_cut=1,inc_score1=1,inc_score2=None,
               U=None,inverse_output=False,rename_columns=False,generateMissing=True):
    self.data=data
    self.input_labels=input_labels
    self.preprocessing=False
    self.n_cut=n_cut
    self.inc_score1=inc_score1
    self.inc_score2=inc_score2
    self.U=U
    self.inverse_output=inverse_output
    self.rename_columns=rename_columns
    self.case_col=case_col
    self.output_labels=output_labels
    self.generate_missing = generateMissing
    self.prime_implicants = None
    self.prepare_rows_called = False
    self.multi_output=len(self.output_labels)>1
    self.rename_dictionary=None
    
  
  def preprocess_data(self):
    data_tmp=self.data.copy()
    if self.input_labels is None:
          self.input_labels=[x for x in list(self.data.columns) if (x not in self.output_labels) and (x!=self.case_col)]
    
    
    if(self.case_col is None or self.case_col == '-None-'):
        data_tmp['case_col']=data_tmp.index.values
        self.case_col = 'case_col'
    if (self.inverse_output):
        data_tmp[self.output_labels]=abs(data_tmp[self.output_lables]-1)
    
    params = {'Inc_{}'.format(c) : (c,inclusion_score) for c in self.output_labels}
    data_grouped=data_tmp.groupby(self.input_labels).agg(n=(self.case_col, 'count'), 
                                                 Cases=(self.case_col, concatenate_strings),
                                                 **params)
    data_grouped=data_grouped[data_grouped['n']>=self.n_cut]
    inc_columns=['Inc_{}'.format(i) for i in self.output_labels]
    if(self.inc_score2 is None):
        data_grouped[self.output_labels]=(data_grouped[inc_columns]>=self.inc_score1).astype(int)
    else:
        if(self.U is None):
            raise Exception('When inc.score2 is specified, U must be specified as well.')
        if(self.U != 0 and self.U != 1):   
            raise Exception('U must be 0 or 1.')
        if (self.U == 1):
            data_grouped[self.output_labels]=(data_grouped[inc_columns]> self.inc_score1).astype(int)
        if (self.U == 0):   
            data_grouped[self.output_labels]=(data_grouped[inc_columns]>= self.inc_score2).astype(int)
            
        
    res = data_grouped.reset_index()
    if self.rename_columns:
        rename_dic = {k:v for k,v in zip(self.input_labels, COLUMN_LABELS[:len(self.input_labels)])}
        self.rename_dictionary=rename_dic
        res.columns = map(lambda x: rename_dic[x] if x in rename_dic.keys() else x, res.columns)
        l=len(self.input_labels)
        self.input_labels=COLUMN_LABELS[:l]

    self.preprocessed_data_raw=res
    self.preprocessed_data=res[self.input_labels+self.output_labels]   
    self.preprocessing=True
   
  def get_rename_dictionary(self):
      return self.rename_dictionary
      
  def get_preprocessed_data(self):
              
      if not self.preprocessing:
          self.preprocess_data()
      return self.preprocessed_data
    
 
  def prepareRows(self):
      if not self.preprocessing:
          self.preprocess_data()
      
      mask1=self.preprocessed_data[self.output_labels]
      nr_rows=len(self.preprocessed_data.index)
      n=len(mask1.columns)
      non_zero_output=[1]*n
      multi_mask=self.preprocessed_data[self.output_labels].isin(non_zero_output)
      mask=multi_mask.aggregate(any, axis=1)
      positiveRows=self.preprocessed_data[mask]
      columns=[col for col in positiveRows.columns if col not in self.output_labels]
      positiveInputs=positiveRows[columns]
      positiveInputs_rownames=list(positiveInputs.index)
      inputs=self.preprocessed_data.drop(self.output_labels,axis=1)
      if len(self.input_labels)==1:
          dim_corrected=[len(set(inputs))]
      else:
          dim=inputs.apply(lambda x: pd.unique(x).tolist(),axis=0, result_type='reduce').array
          
          dim_corrected=[]
          for ar in dim:
              if len(ar)>1:
                  dim_corrected.append(ar)
              else:
                  dim_corrected.append([0,1])
      levels=[len(x) for x in dim_corrected]
      cares_indexes=list()

      if (self.generate_missing):
        allInputs=pd.DataFrame(itertools.product(*dim_corrected))
        zero_output=[0]*n
        multi_mask_zero=self.preprocessed_data[self.output_labels].isin(zero_output)
        mask_zero=multi_mask_zero.aggregate(all,axis=1)
        negativeRows=self.preprocessed_data[mask_zero]
        negativeInputs= negativeRows[columns]
        indexes=list()
        for x in negativeInputs.values:
          ind=find_index2(allInputs, x)
          indexes.append(ind[0])
  
        allInputs_table=allInputs.drop(indexes)
        allInputs_table.columns = columns


        for x in positiveInputs.values:
          ind=find_index2(allInputs_table,x)
          cares_indexes.append(ind[0])
        
        if self.multi_output:
            outcols_merge = pd.merge(allInputs_table, positiveRows, how='inner',on=columns)
            outcols=outcols_merge[self.output_labels]
            self.outputcolumns =outcols.transpose().to_numpy()
        else:
            self.outputcolumns=[1]*nr_rows

        self.levels = levels
        self.cares = cares_indexes
        self.table = allInputs_table.to_numpy()
        #self.outputcolumns =outcols.transpose().to_numpy()
        self.labels = columns
        self.positive_cares = positiveInputs_rownames
      else:
        cares_indexes=[x for x in range(0,len(positiveInputs.to_numpy()))]
        self.levels = levels
        self.cares = cares_indexes
        self.table = positiveInputs.to_numpy()
        self.outputcolumns = positiveRows[self.output_labels].transpose().to_numpy()
        self.labels = columns
        self.positive_cares = positiveInputs_rownames
      self.prepare_rows_called = True
    
  


  def get_prime_implicants(self):
    if self.prime_implicants is not None:
      return self.prime_implicants
  
    if not self.prepare_rows_called:
        self.prepareRows()
        
    if len(self.table) == 0:
        self.prime_implicants = tuple()
        return self.prime_implicants


    table = self.table.astype(int).tolist()
    column_number = len(table[0])
    preprocessed_table=preprocess_input(table)
    prime_implicants=[]
    groups=create_groups(preprocessed_table,column_number, self.cares,self.outputcolumns,self.multi_output)
    reduction_nr = 0
    while(True):
        reduction_res=reduction_step(groups,column_number,self.multi_output)
        if((reduction_res)['implicants']):
            for i in reduction_res['implicants']:
                prime_implicants.append(i)
        if (reduction_res['reduction']is False):
            break
        groups=reduction_res['groups']
        reduction_nr += 1
            
    prime_implicants = eliminate_minterms(table, prime_implicants, self.levels,self.multi_output)
    coverage_dict = create_care_translation_dict(self.cares, self.positive_cares)

    if self.multi_output:
         self.prime_implicants=tuple(Implicant_multi(minterm_to_str(x[0], self.levels, self.labels,0,self.multi_output),{coverage_dict[y] for y in x[1]},outputs=list(x for x in x[2])) for x in prime_implicants)
    else:
         self.prime_implicants =  tuple(Implicant(minterm_to_str(x[0], self.levels, self.labels,0,self.multi_output),x[0],{coverage_dict[y] for y in x[1]}) for x in prime_implicants)
      
    return self.prime_implicants


  def coverage_matrix(self):
    
    if not self.get_prime_implicants():
        self.get_prime_implicants()
    cares=set().union(*[set(x.coverage) for x in self.prime_implicants])
   
    res= np.zeros((len(cares), len(self.prime_implicants)), dtype=bool)
    res_idx_to_care = {v:k for k,v in enumerate(cares)}
    for row_nr,implicant in enumerate(self.prime_implicants):
      for x in implicant.coverage:
        if x not in cares:
          continue
        column_nr = res_idx_to_care[x]
        res[column_nr, row_nr] = True
    if self.multi_output:
      return pd.DataFrame(res.transpose(),columns=cares,index=['{}, {}'.format(x.implicant,x.outputs) for x in self.prime_implicants])    
    return pd.DataFrame(res.transpose(),columns=cares,index=[(x.implicant) for x in self.prime_implicants])
    

  

       
  def get_irredundant_sums(self, max_depth = None):
   if not self.get_prime_implicants():
        self.get_prime_implicants()
        
   if len(self.prime_implicants) == 0:
       return []
        
   if self.multi_output:
        raise RuntimeError("irredudant sums are not supported in multi output mode. Use get_irredundant_systems")
   result=find_irredundant_sums(([(i, i.coverage) for i in self.prime_implicants]),self.cares,max_depth)
   irredundant_objects=[]
   for i,system in enumerate(result):
	   irredundant_objects.append(Irredundant_system(system,i+1))
   return irredundant_objects   
    
  def get_irredundant_systems(self):
   if not self.get_prime_implicants():
        self.get_prime_implicants()
   if not self.multi_output:
        raise RuntimeError("irredudant systems are not supported in single output mode. Use get_irredundant_sums")     
   res,l=self._single_ir_systems_for_multi_output()

   mult_input = [set(frozenset(imp for imp in irs) for irs in f) for f in res]
   reduction_result = reduce(boolean_multiply, mult_input)
   res=[]
   index=0
   for r in reduction_result:
     index=index+1  
     single_res = []
     for j in range(l):
       single_res.append([self.prime_implicants[i] for i in r if j+1 in self.prime_implicants[i].outputs])
     res.append(Irredundant_systems_multi(single_res,index,self.output_labels))
   return res



  def _single_ir_systems_for_multi_output(self):
    l=len(self.output_labels)
    imp_per_output=[]
    for i in range(l):
      imp_per_output.append(list())
    for i,impl in enumerate(self.prime_implicants):
      for j in range(l):
        if (j+1) in set(impl.outputs):
          imp_per_output[j].append((i,impl)) 
    res=[]
    for k,system in enumerate(imp_per_output):
      coverage=set().union(*[set(i[1].coverage) for i in system])
      result=find_irredundant_sums(([(i, impl.coverage) for i,impl in system]),coverage)
      irredundant_objects=[]
      for i,system in enumerate(result):
          irredundant_objects.append(Irredundant_system(system,i+1))
       
      for i in irredundant_objects:
          res.append([i.system for i in irredundant_objects])
 
    return res,l
 


class Irredundant_systems_multi():
  def __init__(self,system_multiple,index,output_labels):
      self.system_multiple=system_multiple
      self.index=index
      self.output_labels=output_labels
  
  def __str__(self):
      res=""
      res+='---- Solution {} ----\n'.format(self.index)
      for j, system in enumerate(self.system_multiple):
         if any(str(impl.implicant) == '1' for impl in system):
             res+=('1 <=> {}\n'.format(self.output_labels[j]))
         elif system==[]:
             res+=('0 <=> {}\n'.format(self.output_labels[j]))
         else:
             res+=('{1} <=> {0}\n'.format(self.output_labels[j], ' + '.join(impl.implicant for impl in system)))
      res+='\n'
      return res
  def __repr__(self):
      return str(self)
      
 

class Irredundant_system():
  def __init__(self,system,index):
	  self.system=system
	  self.index=index
      #self.raw_implicants
    
  def __str__(self):
     return 'M{}:{}'.format(self.index,' + '.join(str(i.implicant) for i in self.system))
 
  def impl_coverag(self):
      res = {}
      for i,impl_i in enumerate(self.system):
          cov_out=set()
          for j,impl_j in enumerate(self.system):
              if j!=i:
                  cov_out.update(impl_j.coverage)
          cov_in=impl_i.coverage-cov_out

          res[str(impl_i.implicant)] = len(cov_in) / len(cov_in.union(cov_out))
      return res
                  

  def __repr__(self):
     return str(self)
 
  def nr_implicants(self):
     return len(self.system)
 
  def coverage_score(self, data, input_columns, output_column):
        tmp_data = data[data[output_column]==1][input_columns]
        return tmp_data.apply(
            lambda row_series: 1.0 if any(all(x in y for x,y in zip(row_series.values, i.raw_implicant)) for i in self.system) else 0.0, axis = 1).mean()
 
  def inclusion_score(self,data, input_columns, output_column):
       tmp_data = data[input_columns]
       mask = tmp_data.apply(
            lambda row_series: any(all(x in y for x,y in zip(row_series[input_columns].values, i.raw_implicant)) for i in self.system), axis = 1)
       return data.loc[mask,output_column].mean()
      



def set_to_str(s,levels,label, is_multi_level):
    if len(s)==levels:
        return ''

    if not is_multi_level:
        if 0 in s:
            return label.lower()
        else:
            return label.upper()
    if len(s)==1:
        return '{}{{{}}}'.format(label,tuple(s)[0])
    return '{}{{{}}}'.format(label,','.join(str(x) for x in s))

def minterm_to_str(minterm, levels, labels, tag,multi_output):
    is_multi_level = any(x > 2 for x in levels)
    tmp = [set_to_str(x, y, z, is_multi_level) for x,y,z in zip(minterm, levels, labels)]
    res = '{}'.format('*'.join(x for x in tmp if x != ''))
    return res if res != '' else '1'



class Multiple_output_item:
     
    def __init__(self, minterm, coverage, tag):
        self.minterm = tuple(x for x in minterm)
        self.coverage = frozenset(x for x in coverage)
        self.is_reduced=False
        self.tag=frozenset(x for x in tag)
        
    def __str__(self):
        return('{0}, {1},{2}, tag={3}'.format(
                                    self.minterm, 
                                    self.coverage,
                                    self.is_reduced,
                                    self.tag))
    def __repr__(self):
        return str(self)
        
    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (self.minterm == other.minterm and 
                self.coverage == other.coverage and 
                self.is_reduced == other.is_reduced and
                self.tag == other.tag)
                
    def __hash__(self):
        return hash((self.minterm, self.coverage, self.is_reduced))
    
    def can_be_reduced(self,other):
        if len(self.tag.intersection(other.tag)) == 0:
            return False
        n=len(self.minterm)
        diff=0
        for i in range(0,n):
            if self.minterm[i]!=other.minterm[i]:
                diff=diff+1
        if diff==1:
            return(True)
        else:
            return(False)
    
    def reduce(self, other):
        if not self.can_be_reduced(other):
            return None
        new_tag=self.tag.intersection(other.tag)
        if new_tag ==self.tag:
            self.is_reduced=True
        if new_tag==other.tag:
            other.is_reduced = True
        n=len(self.minterm)
        new_minterm=list()
        for i in range(0,n):
            new_minterm.append({})
        for i in range(0,n):
            if(self.minterm[i]==other.minterm[i]):
                new_minterm[i]=self.minterm[i]
            else:
                new_minterm[i]=self.minterm[i].union(other.minterm[i])
        return Multiple_output_item(new_minterm, 
                       self.coverage.union(other.coverage),new_tag)


class Implicant:

    def __init__(self,implicant,raw_implicant,coverage,cov_score=None,cov_u=None,incl_score=None):
      self.implicant=implicant
      self.raw_implicant = raw_implicant
      self.coverage=coverage
      
    def __str__(self):
        return('{0}:{1}'.format(self.implicant, self.coverage))

    def __repr__(self):
        return str(self)
    
    def coverage_score(self, data, input_columns, output_column):
        if (len(input_columns) != len(self.raw_implicant)):
            raise RuntimeError(
                'Size of input columns ({}) does not match implicant size({})'.format(len(input_columns), 
                                                                                      len(self.raw_implicant)))
        tmp_data = data[data[output_column]==1][input_columns]
        return tmp_data.apply(
            lambda row_series: 1.0 if all(x in y for x,y in zip(row_series.values, self.raw_implicant)) else 0.0, axis = 1).mean()
    
    def inclusion_score(self, data, input_columns, output_column):
        if (len(input_columns) != len(self.raw_implicant)):
            raise RuntimeError(
                'Size of input columns ({}) does not match implicant size({})'.format(len(input_columns), 
                                                                                      len(self.raw_implicant)))
        tmp_data = data[input_columns]
        return tmp_data.apply(
            lambda row_series: 1.0 if all(x in y for x,y in zip(row_series.values, self.raw_implicant)) else 0.0, axis = 1).sum() / len(tmp_data.index)
        
        

class Implicant_multi:

    def __init__(self,implicant,coverage,outputs,cov_score=None,cov_u=None,incl_score=None):
      self.implicant=implicant
      self.coverage=coverage

      self.outputs=outputs

    def __str__(self):
        return('{0}:{1},{2}'.format(self.implicant, self.coverage, self.outputs))

    def __repr__(self):
        return str(self)


class Multi_value_item:
    
    def __init__(self, minterm, coverage):
        self.minterm = tuple(x for x in minterm)
        self.coverage = frozenset(x for x in coverage)
        self.is_reduced=False
        
    def __str__(self):
        return('{0}, {1},{2}'.format(
                                    self.minterm, 
                                    self.coverage,
                                    self.is_reduced
                                    ))
    def __repr__(self):
        return str(self)
        
    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (self.minterm == other.minterm and 
                self.coverage == other.coverage and 
                self.is_reduced == other.is_reduced)
                
    def __hash__(self):
        return hash((self.minterm, self.coverage, self.is_reduced))
        
    def can_be_reduced(self,other):
        n=len(self.minterm)
        diff=0
        for i in range(0,n):
            if self.minterm[i]!=other.minterm[i]:
                diff=diff+1
        if diff==1:
            return(True)
        else:
            return(False)
    
    def reduce(self, other):
        if not self.can_be_reduced(other):
            return None
        self.is_reduced=True
        other.is_reduced=True
        n=len(self.minterm)
        new_minterm=list()
        for i in range(0,n):
            new_minterm.append({})
        for i in range(0,n):
            if(self.minterm[i]==other.minterm[i]):
                new_minterm[i]=self.minterm[i]
            else:
                new_minterm[i]=self.minterm[i].union(other.minterm[i])
        return Multi_value_item(new_minterm, 
                       self.coverage.union(other.coverage))














