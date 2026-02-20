import cora
import pandas as pd

# binary data example
alg_to_use = "ON-DC"
# uncomment line below to use ON-OFF algorithm
# alg_to_use = 'ON-OFF'

data1 = pd.read_csv("GrossCarvin2011.csv", sep=";")

# create OptimizationContext
context1 = cora.OptimizationContext(
    data1, output_labels=["TORT"], case_col="Case", algorithm=alg_to_use
)


# calculate the prime implicants with ON-DC algorithm
context1.get_prime_implicants()

# get the prime implicants details
context1.pi_details()

# get prime implicant chart
context1.prime_implicant_chart()

# calculate irredundant sums
context1.get_irredundant_sums()

# get the summary table
context1.get_solution_dataframe()
