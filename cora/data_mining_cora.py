import itertools
import pandas as pd
import cora



class TupleResult:
    
    def __init__(self, combination,irrendudant_sums, score):
        
        self.combination=combination
        if(len(irrendudant_sums) == 1 and (str(
                irrendudant_sums[0].system[0].implicant)=='1')):
            self.nr_irr_sums = 0
            self.min_nr_pi = 0
            self.max_nr_pi = 0
            self.score=score = 0
        else:     
            self.nr_irr_sums = len(irrendudant_sums)
            self.min_nr_pi = (min(x.nr_implicants() for x in irrendudant_sums)
                               if self.nr_irr_sums > 0 else 0)
            self.max_nr_pi = (max(x.nr_implicants() for x in irrendudant_sums)
                               if self.nr_irr_sums > 0 else 0)
            self.score=score
 
        
    

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

        ir_sums = data_object.get_irredundant_sums()
        

        if len(ir_sums) > 0:
          score = round(max(x.inclusion_score()*x.coverage_score() for x
                                                  in ir_sums),3)
          tr = TupleResult([x for x in comb],ir_sums, score)
          res.append(tr)
    
        else:
          res.append(TupleResult(comb, ir_sums, 0))
    
    rows = [(x.combination,
                     x.nr_irr_sums,
                     x.max_nr_pi, 
                     x.min_nr_pi,
                     x.score) for x in res]
    result = pd.DataFrame(rows,
                          columns = ['Combination',
                                     'Nr_of_sums',
                                     'Max_PIs',
                                     'Min_PIs',
                                     'Score'])
      
    return result


    

