import itertools
import pandas as pd
import cora



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

"""

Parameters

----------

data : dataframe

output_labels : an array of strings
    The names of the outcome columns from the data frame.
    
len_of_tuple : int
    Number indicating how many variables from the original data 
    are used in the computation.

case_col : string
    The name of the column from the data frame containing the case ids.

n_cut : int
    The minimum number of cases under which a truth table row is declared as a
    remainder
        
inc_score1 : float
	The minimum sufficiency inclusion score for an output function value of "1"
    
inc_score2 : float
	The maximum sufficiency inclusion score for an output function value of "0"

U : int
    The U number is either 0 or 1.

algorithm : string
            ON-DC or ON-OFF
        

Returns
-------

result : dataframe
    
    
"""   

def data_mining(data,
                output_labels,
                len_of_tuple,
                case_col=None,
                n_cut=1,
                inc_score1=1,
                inc_score2=None,
                Uvalue=None,
                algorithm="ON-DC"):
    res = []

        
    for i,comb in enumerate(itertools.combinations([x for x in data.columns
   
                                                if (x not in output_labels and
                                                        x!=case_col)],
                                                       len_of_tuple)):
        cols = list(comb)
        data_object = cora.OptimizationContext(data,
                                               output_labels,
                                                cols,
                                                case_col=case_col,
                                                n_cut=n_cut,
                                                inc_score1=inc_score1,
                                                inc_score2=inc_score2,
                                                U=Uvalue,
                                                algorithm=algorithm)
            
        if len(output_labels) == 1:
            ir_sys = data_object.get_irredundant_sums()
        else:
            ir_sys = data_object.get_irredundant_systems()
            
    
        if len(ir_sys) > 0:
                inc_score = round(max(x.inclusion_score() for x in ir_sys),3)
                cov_score = round(max(x.coverage_score() for x in ir_sys),3)
                score = round(max(x.inclusion_score()*x.coverage_score()
                                                           for x  in ir_sys),3)
                tr = TupleResult(cols,
                                 len(output_labels),
                                 ir_sys,
                                 inc_score,
                                 cov_score,
                                 score)
                res.append(tr)
        
        else:
                res.append(TupleResult(cols, len(output_labels), ir_sys, 0, 0, 0))
            
        
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


if __name__ == '__main__':
    df = pd.DataFrame([[1,1,0,1,1,1],[0,1,1,1,0,1],[0,0,1,0,0,1],[1,0,1,1,0,1]],
                      columns = ["A","B","C","D","Z","P"])
    print(data_mining(df,['Z','P'],2))