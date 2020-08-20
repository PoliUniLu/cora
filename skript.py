#!/usr/bin/env python3
# -*- coding: utf-8 -*-s

import pandas
import cora
#import itertools

data=pandas.read_csv("//Users/zuzka/Documents/git/SwissMinaret.csv",sep=';')
#data=pandas.read_csv("//Users/zuzka/Documents/git/new_samples/kluver_agg.csv",sep=',')

#tmp=cora.preprocess_data(data,out_cols=["M"])
#data_new=tmp[["A","L","S","T","X","M"]]
#data_new=data_new.sort_values
#print("new_data:{}".format((data_new)))
#print(tmp)
c=cora.Chart(data,["M"],case_col="Case")
k=c.get_prime_implicants()
print('prie_implicants:{}'.format(k))
#print(k[0].coverage_score(data, ["A","L","S","T","X"], "M"))
#print(k[0].inclusion_score(data, ["A","L","S","T","X"], "M"))
#P=c.prime_implicants
z=c.get_irredundant_sums()
print("systems:{}".format(z))
