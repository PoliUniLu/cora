import pandas as pd
import itertools
import cora
data=pd.read_csv("//Users/zuzka/Documents/git/new_samples/kluver_agg.csv",sep=',')
print(data)


class TupleResult:
    
    def __init__(self, irrendudant_sums, score):
        self.nr_irr_sums = len(irrendudant_sums)
        self.min_ir_sum = min(x.nr_implicants() for x in irrendudant_sums)
        self.max_ir_sum = max(x.nr_implicants() for x in irrendudant_sums)
        self.score=score
    

def data_mining(data,out_col,len_of_tupple):
    res=[]
    
    for i,comb in enumerate(itertools.combinations([x for x in data.columns if x != out_col], len_of_tupple)):
        cols = list(comb)+[out_col]
        data_tmp=cora.preprocess_data(data[cols],out_col)
        data_chart=cora.Chart(data_tmp[cols],[out_col])
        ir_sums = data_chart.get_irredundant_sums()
        score = max(x.inclusion_score(data, list(comb), out_col)*x.coverage_score(data, list(comb), out_col) for x in ir_sums)
        tr = TupleResult(ir_sums, score)
        res.append(tr)
        print('{} -> nr_sums: {} min nr of Pis: {} max nr of PIs:{} score: {}'.format(comb, 
                                                         tr.nr_irr_sums, 
                                                         tr.min_ir_sum,
                                                         tr.max_ir_sum,
                                                         score
                                                        )) 
    return res

data_mining(data, "E", 3)