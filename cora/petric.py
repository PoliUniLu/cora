from functools import reduce
from native_petric import petric
import sys


def _find_irredundant_sums_native(implicants_with_coverage, coverage):
    
    petric_input  = [x[1] for x in implicants_with_coverage]
    #print('Petric called with {} immplicants'.format(len(petric_input)))
    #print(petric_input)
    #sys.exit(0)
    _ , sums = petric(petric_input)
   # print('Petic done, foudn {} sums'.format(len(sums)))
    
    res = []
    for x in sums:
        res.append(implicants_with_coverage[i][0] for i in x)
    return res
    
    

def _find_irredundant_sums(implicants_with_coverage, coverage, max_depth = None):
	if max_depth is None:
		max_depth = len(implicants_with_coverage)

	results = []
	_find_irrendundant_sums_internal(implicants_with_coverage, [],
									 set(), len(coverage), results,
									 max_depth)
	return results

def _find_irrendundant_sums_internal(implicants_with_coverage, partial_solution,
									 coverage, to_cover, all_solutions, max_depth):

	# If we reached maximal depth / maximal length of the sum, do not continue.
	if len(partial_solution) > max_depth:
		return

	if not _is_irredundant_sum(partial_solution, coverage):
			return

	# If the coverage is empty, this means the partial solution already covers everything.
	if len(coverage) >= to_cover:
		all_solutions.append([x[0] for x in partial_solution])
		return

	for i, implicant_with_coverage in enumerate(implicants_with_coverage):
		imp, imp_coverage = implicant_with_coverage
		if len(imp_coverage.difference(coverage)) == 0:
			continue
		new_partial_solution = partial_solution + [implicant_with_coverage]
		new_implicants = [implicants_with_coverage[j] 
		                  for j in range(len(implicants_with_coverage)) 
		                  if j > i]
		new_coverage = coverage.union(imp_coverage)
		_find_irrendundant_sums_internal(new_implicants,
										 new_partial_solution,
										 new_coverage,
										 to_cover,
										 all_solutions,
										 max_depth)

def _is_irredundant_sum(partial_solution_with_coverage, coverage):
	coverage_counts = {k:0 for k in coverage}
	for imp, imp_coverage in partial_solution_with_coverage:
		for x in imp_coverage:
			coverage_counts[x] += 1
	for _, imp_coverage in partial_solution_with_coverage:
		if all(coverage_counts[x] > 1 for x in imp_coverage):
			return False
	return True



def _implicants_coverage(pi_chart):
  final_res=[]
  coverage=set()
  for i,row in enumerate(pi_chart.index):
    coverage = set(index for index, value in pi_chart.iloc[i,:].items() if value)
    final_res.append((row,coverage))
  return final_res, pi_chart.columns
  







def _boolean_multiply(x, y):
    assert(isinstance(x,set))
    assert(isinstance(y,set))
    res = set()
    for x_i in x:
        for y_i in y:
            tmp = x_i.union(y_i)
            # X+XY=X
            # if there is superset of tmp in result,
            # we want to remove it.
            res = set(filter(lambda z: not tmp.issubset(z), res))
            # if tmp is superset of existing term in res,
            # we don't add tmp.
            if all(not z.issubset(tmp) for z in res):
                res.add(tmp)
    return res









