#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 23:26:45 2020

@author: zuzka
"""
import pandas as pd
import cora

data=pd.read_csv("//Users/zuzka/Downloads/T1-T3.csv",sep=",")
print(data)
c=cora.Chart(data,['T1','T3'])
s=c.get_irredundant_systems()
print(s)