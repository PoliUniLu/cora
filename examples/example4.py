import cora
import pandas as pd

#data mining example
alg_to_use = 'ON-DC'

#uncomment line below to use ON-OFF algorithm
#alg_to_use = 'ON-OFF'

data1 = pd.read_csv("McCluskeyF1F2.csv",sep=';')

cora.data_mining(data1,
                len_of_tuple=2,
                output_labels=["F1","F2"],
                algorithm = alg_to_use)



