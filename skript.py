#!/usr/bin/env python3
# -*- coding: utf-8 -*-s

import pandas
import cora
#import itertools

data=pandas.read_csv("//Users/zuzka/Documents/cora/test_multi.csv",sep=';')
#data=pandas.read_csv("//Users/zuzka/Documents/git/new_samples/kluver_agg.csv",sep=',')
print(data)

c=cora.Chart(data,["W","X","Y","Z"])
k=c.get_prime_implicants()
#print('prie_implicants:{}'.format(k))
#print(k[0].coverage_score(data, ["A","L","S","T","X"], "M"))
#print(k[0].inclusion_score(data, ["A","L","S","T","X"], "M"))
#P=c.prime_implicants
z=c.get_irredundant_systems()
print("systems:{}".format(z))
