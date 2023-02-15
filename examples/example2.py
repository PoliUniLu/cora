import cora
import pandas as pd

#multiple outputs example

alg_to_use = 'ON-OFF'
#uncomment line below to use ON-OFF algorithm
#alg_to_use = 'ON-DC'

data1 = pd.read_csv("SwissMinaret.csv",sep=',')
# create OptimizationContext
context1 = cora.OptimizationContext(data1,
                                   output_labels=["X","M"],
                                   algorithm = alg_to_use)



# compose the truth table
context1.get_preprocessed_data()

#calculate the prime implicants
context1.get_prime_implicants()

# calculate irredundant systems
print(context1.get_irredundant_systems())

# get the summary table
context1.get_solution_dataframe()




