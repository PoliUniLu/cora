import itertools
import pandas as pd
import cora
import re

TEMPORAL_COL_PATTERN = re.compile("^([a-zA-Z0-9-\_]+)\{([0-9]+(,[0-9]+)*)\}$")                                              


class TupleResult:
    
    def __init__(self,
                 combination,
                 out_len,
                 irrendudant_systems,
                 inc_score,
                 cov_score,
                 score):
        
        self.combination=combination
        if(out_len ==1 and len(irrendudant_systems) == 1 and (str(
                irrendudant_systems[0].system[0].implicant)=='1')):
            self.nr_irr_systems = 0
            self.inc_score = 0
            self.cov_score = 0
            self.score=score = 0
        else:     
            self.nr_irr_systems = len(irrendudant_systems)
            self.inc_score = inc_score
            self.cov_score = cov_score
            self.score = score

#valid_comb: ((1, 'PUBLIC'), (2, 'ELITE')) , ['E_BEFORE_A{2,3}']

def valid_combination(combination,temp_cols):
    #combination_cols  = set(x[1] for x in combination) 
    #temp_cols_string = {m.group(1) 
     #                   for m in [TEMPORAL_COL_PATTERN.match(i)
     #                   for i in temp_cols]}
    if len(temp_cols)<1: 
       return True
 
    else:
        cols_indexes = set(x[0] for x in combination)
        tmp_t_labels = {m.group(2) 
                        for m in [TEMPORAL_COL_PATTERN.match(i)
                        for i in temp_cols]}
        temp_index_tmp = [x.split(",") for x in tmp_t_labels]
        temp_index = {int(x) for y in temp_index_tmp for x in y}
        
        return temp_index.issubset(cols_indexes)
 

    
    
"""

Parameters

----------

data : dataframe

out_col : an array of strings
    The names of the outcome columns from the data frame.
    
len_of_tupple : int
    Number indicating how many variables from the original data 
    are used in the computation.

case_col : string
    The name of the column from the data frame containing the case ids.

cut : int
    The minimum number of cases under which a truth table row is declared as a
    remainder
        
inc1 : float
	The minimum sufficiency inclusion score for an output function value of "1"
    
inc2 : float
	The maximum sufficiency inclusion score for an output function value of "0"

U : int
    The U number is either 0 or 1.


Returns
-------

result : dataframe
    
    
"""   

def data_mining(data,
                out_col,
                len_of_tupple,
                case_col=None,
                temp_cols=None,
                cut=1,
                inc1=1,
                inc2=None,
                Uvalue=None):
    res = []
    
    for i,comb in enumerate(itertools.combinations([(ind+1,x) for ind,x in 
                                                    enumerate(data.columns) 
                                                    if (x not in out_col and
                                                    
                                                        x!=case_col)
                                                    ], 
                                                   len_of_tupple)):
       
        cols = list(x[1] for x in comb)

        tmp_temp = [x for x in temp_cols if any(x.startswith(y+'{') for y in cols)]
        input_cols = [x for x in cols if all(not y.startswith(x+'{') for y in tmp_temp)]
        if valid_combination(comb,tmp_temp):
            
            
            if len(tmp_temp) >= 1 :
                temp_labels = tmp_temp
            else:
                temp_labels = None
      
            
            
            data_object = cora.OptimizationContext(data,
                                                   out_col,
                                                   input_cols,
                                                   case_col=case_col,
                                                   temporal_labels= temp_labels,
                                                   n_cut=cut,
                                                   inc_score1=inc1,
                                                   inc_score2=inc2,
                                                   U=Uvalue)
            if len(out_col) == 1:
                ir_sys = data_object.get_irredundant_sums()
            else:
                ir_sys = data_object.get_irredundant_systems()
            
    
            if len(ir_sys) > 0:
                inc_score = round(max(x.inclusion_score() for x in ir_sys),3)
                cov_score = round(max(x.coverage_score() for x in ir_sys),3)
                score = round(max(x.inclusion_score()*x.coverage_score() 
                                                           for x  in ir_sys),3)
                tr = TupleResult(cols,
                                     len(out_col),
                                     ir_sys,
                                     inc_score,
                                     cov_score,
                                     score)
                res.append(tr)
        
            else:
                    res.append(TupleResult(cols, len(out_col), ir_sys, 0, 0, 0))
            
        
    rows = [(x.combination,
                 x.nr_irr_systems,
                 x.inc_score,
                 x.cov_score,
                 x.score) for x in res]
    result = pd.DataFrame(rows,
                              columns = ['Combination',
                                         'Nr_of_systems',
                                         'Inc_score',
                                         'Cov_score',
                                         'Score'])
          
    return result


    
