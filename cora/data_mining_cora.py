import itertools
import cora



class TupleResult:
    
    def __init__(self, combination,irrendudant_sums, score):
        self.combination=combination
        self.nr_irr_sums = len(irrendudant_sums)
        self.min_ir_sum = min(x.nr_implicants() for x in irrendudant_sums) if self.nr_irr_sums > 0 else 0
        self.max_ir_sum = max(x.nr_implicants() for x in irrendudant_sums) if self.nr_irr_sums > 0 else 0
        self.score=score
    
    

def data_mining(data,out_col,len_of_tupple,case_col=None,cut=1,inc1=1,inc2=None,Uvalue=None):
    res=[]
    
    for i,comb in enumerate(itertools.combinations([x for x in data.columns if (x not in out_col and x!=case_col)], len_of_tupple)):
        cols = list(comb)
        data_object=cora.Chart(data,out_col,cols,case_col=case_col,n_cut=cut,inc_score1=inc1,inc_score2=inc2,U=Uvalue)


        ir_sums = data_object.get_irredundant_sums()
        #print(data_object.prime_implicants)
        if len(ir_sums) > 0:
          score = round(max(x.inclusion_score(data, list(comb), out_col[0])*x.coverage_score(data, list(comb), out_col[0]) for x in ir_sums),3)
          tr = TupleResult(comb,ir_sums, score)
          res.append(tr)
        else:
          res.append(TupleResult(comb, ir_sums, 0))
      
    return res


    

