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
   

def data_mining(data,
                out_col,
                len_of_tupple,
                case_col=None,
                cut=1,
                inc1=1,
                inc2=None,
                Uvalue=None):
    res = []
    
    for i,comb in enumerate(itertools.combinations([x for x in data.columns 
                                                    if (x not in out_col and
                                                        x!=case_col)], 
                                                   len_of_tupple)):
        cols = list(comb)
        data_object = cora.OptimizationContext(data,
                                               out_col,
                                               cols,
                                               case_col=case_col,
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
            tr = TupleResult([x for x in comb],
                                 len(out_col),
                                 ir_sys,
                                 inc_score,
                                 cov_score,
                                 score)
            res.append(tr)
    
        else:
                res.append(TupleResult(comb, len(out_col), ir_sys, 0, 0, 0))
        
    
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


    

