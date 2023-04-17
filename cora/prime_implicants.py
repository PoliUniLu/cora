"""
This module contains the implementation of Boolean minimization algorithms.
"""

import numpy as np
import pandas as pd
import itertools
from functools import reduce
import string
import re
from collections import defaultdict

from .petric import _find_irredundant_sums_native, _boolean_multiply
from .on_off_alg import _on_off_grouping, _reduction
from .on_off_mo_alg import on_off_grouping_mo, reduction_mo
from .multiply import _transform_to_raw_implicant
from .draft_cubes import _find_implicants_cubes, _transform_to_raw_imp
from .essential import get_essential_implicants, _transform_to_raw_impl, \
    _reduce_the_onset
COLUMN_LABELS = list(string.ascii_uppercase) + ["AA", "BB", "CC", "DD",
                                                "EE", "FF"]
OUTPUT_PATTERN = re.compile("^([a-zA-Z0-9]+)\{([0-9]+(,[0-9]+)*)\}$")
REGULAR_OUTPUT = re.compile("^([a-zA-Z0-9]+)")


def _concatenate_strings(arr):
    return ','.join([str(x) for x in arr])


def _inclusion_score(arr):
    return round(sum(arr) / len(arr), 2)


def _count_non_zeros(minterm):
    n = len(minterm)
    non_zeros = 0
    for i in range(0, n):
        if minterm[i] == {0}:
            continue
        else:
            non_zeros = non_zeros + 1
    return non_zeros


def _vector_to_minterm(v):
    return tuple(frozenset({x}) for x in v)


def _preprocess_input(table):
    return [_vector_to_minterm(x) for x in table]


def _create_care_translation_dict(cares, labels):
    return {k: v for k, v in zip(cares, labels)}


def _check_element_coverage(row, element):
    for r, e in zip(row, element[0]):
        if r not in e:
            return False
    return True


def _is_minterm_subset(m1, m2):
    for x, y in zip(m1[0], m2[0]):
        if not x.issubset(y):
            return False
    return True


# Funtion which groups minterms according to the number of non-zero digets

# Paramters: table - the truth table
#            column_number - defines the number of possible non-zero digits
#            cares - indexes of the positive mintermsin in the data
#            output_columns - binary output columns from the data
#            multi_output - a boolean variable
# Output: (a) an array res of objects from classes Multi_value_item or
#             Multiple_output_minterm
#         (b) the i-th element of an array consisting of
#              objects whose minterms have exactly i non-zero digits

# Definition: a minterm of n variables is a product of the variables
# 	      in which each one appears exactly once in true or complemented form.

def _create_groups(table, column_number, cares, outputcolumns, multi_output):
    res = []
    for i in range(0, column_number + 1):
        res.append([])
    index = 0
    if not multi_output:
        for row in table:
            non_zeros = _count_non_zeros(row)
            res[non_zeros].append(MultiValueMinterm(row, set([index])
            if index in cares else set()))
            index = index + 1
        return res
    else:
        nr_of_outputs = len(outputcolumns)
        care_index = 0
        for row in table:
            non_zeros = _count_non_zeros(row)
            if index in cares:
                tag = {i for i in range(1, nr_of_outputs + 1)
                       if outputcolumns[i - 1][care_index] == 1}
                care_index = care_index + 1
            else:
                tag = {i for i in range(1, nr_of_outputs + 1)}
            if len(tag) == 0:
                index = index + 1
                continue
            res[non_zeros].append(MultipleOutputMinterm(row, set([index])
            if index in cares else set(), tag))
            index = index + 1
        return res


# Function performing 1 step elimination of 2 minterms

# Function eliminates minterms first in adjacent groups, second in
# between the group. In both cases minterms differ in only one sigle digit.
# In addition, the function collects implicants, which can't be further
# eliminated to an array- the prime implicants.

# Parameters: groups - an array returned by fun create_groups
#             n - the number of rows in the truth table
#             multi_output - a boolean variable

# Output: dictionary, which contains
# 	  (a) new groups of minterms after the elimination
#         (b) implicants
#         (c) information if any reduction was performed in this single
#             elimination

def _reduction_step(groups, n, multi_output):
    was_any_reduction = False
    new_groups = []
    for i in range(n + 1):
        new_groups.append(set())

    for i in range(0, n):
        if len(groups[i]) == 0:
            continue
        for element1 in groups[i]:
            for element2 in groups[i + 1]:
                if element1.can_be_reduced(element2):
                    reduced_element = element1.reduce(element2)
                    new_groups[i].add(reduced_element)
                    was_any_reduction = True

    for i in range(0, n + 1):
        for element1 in groups[i]:
            for element2 in groups[i]:
                if element1.can_be_reduced(element2):
                    reduced_element = element1.reduce(element2)
                    new_groups[i].add(reduced_element)
                    was_any_reduction = True

    final_implicants = set([])
    for i in range(0, n + 1):
        for element1 in groups[i]:
            if (element1.is_reduced is False and (len(element1.coverage) > 0)):
                if multi_output:
                    final_implicants.add((
                        tuple(element1.minterm),
                        frozenset(element1.coverage),
                        frozenset(element1.tag)))
                else:
                    final_implicants.add((tuple(element1.minterm),
                                          frozenset(element1.coverage)))

    return ({'groups': new_groups,
             'implicants': final_implicants,
             'reduction': was_any_reduction})


# Next set of 3 functions perform an extra elimination step


def _extend_with_first(first, arr, multi_output):
    if multi_output:
        return [(tuple([first]) + x[0], x[1], x[2]) for x in arr]
    return [(tuple([first]) + x[0], x[1]) for x in arr]


def _decomposition(element, levels, multi_output):
    if len(element[0]) == 0:
        return [element]
    if multi_output:
        tmp = _decomposition((element[0][1:], element[1], element[2]),
                            levels[1:], multi_output)
    else:
        tmp = _decomposition((element[0][1:], element[1]),
                            levels[1:], multi_output)

    if len(element[0][0]) == 1 or len(element[0][0]) == levels[0]:
        return _extend_with_first(element[0][0], tmp, multi_output)
    else:
        res = []
        for x in element[0][0]:
            res.extend(_extend_with_first({x}, tmp, multi_output))
        return res


def _fix_coverage_after_decomposition(table, elements, multi_output):
    new_elements = []
    for e in elements:
        new_coverage = frozenset([x for x in e[1] if
                                  _check_element_coverage(table[x], e)])

        if len(new_coverage) > 0:
            if multi_output:
                new_elements.append(tuple([e[0], new_coverage, e[2]]))
            else:
                new_elements.append(tuple([e[0], new_coverage]))

    return new_elements


def _eliminate_minterms(table, elements, levels, multi_output):
    decomposed = []
    for x in elements:
        decomposed.extend(_decomposition(x, levels, multi_output))

    decomposed = _fix_coverage_after_decomposition(table,
                                                   decomposed,
                                                   multi_output)
    n = len(decomposed)
    was_eliminated = [False] * n
    if multi_output:
        for i in range(n):
            for j in range(i + 1, n):
                if _is_minterm_subset(decomposed[i], decomposed[j]) \
                        and (decomposed[i][2]).issubset(decomposed[j][2]):
                    was_eliminated[i] = True
                elif _is_minterm_subset(decomposed[j], decomposed[i]) \
                        and (decomposed[j][2]).issubset(decomposed[i][2]):
                    was_eliminated[j] = True
        return [x for we, x in zip(was_eliminated, decomposed) if not we]

    for i in range(n):
        for j in range(i + 1, n):

            if _is_minterm_subset(decomposed[i], decomposed[j]):
                was_eliminated[i] = True
            elif _is_minterm_subset(decomposed[j], decomposed[i]):
                was_eliminated[j] = True
    return [x for we, x in zip(was_eliminated, decomposed) if not we]



def _find_index2(arr, x):
    return np.where((arr == x).all(axis=1))[0]


def _set_to_str(s, levels, label, is_multi_level):
    if len(s) == levels:
        return ''

    if not is_multi_level:
        if 0 in s:
            return label.lower()
        else:
            return label.upper()
    if len(s) == 1:
        return '{}{{{}}}'.format(label, tuple(s)[0])
    return '{}{{{}}}'.format(label, ','.join(str(x) for x in s))


def _minterm_to_str(minterm, levels, labels, tag, multi_output):
    is_multi_level = any(x > 2 for x in levels)
    tmp = [_set_to_str(x, y, z, is_multi_level) for x, y, z in
           zip(minterm, levels, labels)]
    res = '{}'.format('*'.join(x for x in tmp if x != ''))

    return res if res != '' else '1'


def _calculate_essential_indexes(prime_implicants):
    essentials = []

    essentials_histogram = defaultdict(lambda: 0)
    for x in prime_implicants:
        for y in x[1]:
            essentials_histogram[y] += 1

    for e in essentials_histogram.keys():
        if essentials_histogram[e] == 1:
            essentials.append(e)
    return set(essentials)


class OptimizationContext:
    """
    Optimization Context contains the data and information to preform the
    computation.
    Parameters
    ----------
    data : dataframe
    output_labels : an array of strings
        The names of the outcome columns from the data frame. If the columns
        contain values requiring to map to a boolean range, the set
        of the values can be added at the end of the name string.
    input_labels : array of strings
        The names of the input columns from the data frame
    case_col : string
        The name of the column from the data frame containing the case ids
    n_cut : int
        The minimum number of cases under which a truth table row is declared as a
        don't care.
    inc_score1 : float
        The minimum sufficiency inclusion score for an output function value of "1".
    inc_score2 : float
        The maximum sufficiency inclusion score for an output function value of "0".
    U : int
        The U number is either 0 or 1.
    rename_columns : boolean
        If true, the column names are renamed as single letters in alphabetic
        order.
    algorithm : string
        The name of the optimization algorithm.
        (a): ON-DC, the classical Quine-McCluskey algorithm operating on
                    positive (ON) and don't care (DC) terms.
        (b): ON-OFF (default), the modified version of McCluskey's algorithm operating on
                    positive and negative (OFF) terms.
    """

    def __init__(self,
                 data,
                 output_labels,
                 input_labels=None,
                 case_col=None,
                 n_cut=1,
                 inc_score1=1,
                 inc_score2=None,
                 U=None,
                 rename_columns=False,
                 algorithm="ON-DC"
                 ):
        self.data = data.copy()
        self.input_labels = input_labels
        self.preprocessing = False
        self.validation = False
        self.n_cut = n_cut
        self.inc_score1 = inc_score1
        self.inc_score2 = inc_score2
        self.U = U
        self.rename_columns = rename_columns
        self.case_col = case_col
        self.output_labels = output_labels
        self.prime_implicants = None
        self.irredundat_sums = None
        self.irredundant_systems = None
        self.prepare_rows_called = False
        self.multi_output = len(self.output_labels) > 1
        self.rename_dictionary = None
        self.multivalue_output = False
        self.details = None
        self.sol_details = None
        self.pi_chart = None
        self.input_data = None
        self.algorithm = algorithm
        self.solution_dataframe = None

    # Function to clean and aggregate data. Removes duplicities and
    # inconsistencies.

    def _data_validation(self):

        class InvalidDataException(Exception):
            pass

        if self.input_labels is None:

                inputs = list(x for x in list(self.data.columns)
                              if x != self.case_col)
        else:
                inputs = self.input_labels
        if (not all(self.data[inputs].apply(
                    lambda row_series: all(isinstance(x, int)
                                           for x in row_series), axis=0))):
                raise InvalidDataException("Invalid data input!")

        # outputs

        if all(OUTPUT_PATTERN.match(i) is not None for i in self.output_labels):
            self.multivalue_output = True

        elif not all(REGULAR_OUTPUT.match(i) is not None for i
                     in self.output_labels):
            raise RuntimeError("Unsupported output entered!")

        if self.multivalue_output:
            outcols_value_map = dict()
            for i in self.output_labels:
                values = set(int(x) for x in
                             OUTPUT_PATTERN.match(i).group(2).split(","))
                outcols_value_map[OUTPUT_PATTERN.match(i).group(1)] = values

            for k in outcols_value_map.keys():
                self.data[k] = self.data[k].apply(lambda x: 1 if x in
                                                            outcols_value_map[
                                                                     k] else 0)
                self.output_labels = [str(k) for k in outcols_value_map.keys()]

            self.output_labels_final = [str(k) + str(set(outcols_value_map[k]))
                                        for k in outcols_value_map.keys()]
        else:
            output_values = pd.unique(self.data[self.output_labels].values.ravel('K'))
            if not np.isin(output_values, [0,1]).all():
                raise RuntimeError("Unsupported output entered!"
                                   "\Please specify the analysing"
                                   "\ output values. )")



            self.output_labels_final = self.output_labels

        if self.input_labels is None:
            self.input_labels = [x for x in list(self.data.columns) if
                                 (x not in self.output_labels)
                                 and (x != self.case_col)]

        input_data = self.data[self.input_labels]

        if self.input_data is None:
            self.input_data = input_data

        if (any(input_data.apply(lambda row_series: True if
        len(row_series.unique()) == 1 else False, axis=0))):
            raise InvalidDataException("Please respecify your input data"
                                       + " constants are not allowed!")

        self.validation = True

    def _preprocess_data(self):

        if not self.validation:
            self._data_validation()

        data_tmp = self.data.copy()
        if (self.case_col is None or self.case_col == '-None-'):
            data_tmp['case_col'] = data_tmp.index.values
            self.case_col = 'case_col'
        # inclusion_score is a function for the data aggregation
        params = {'Inc_{}'.format(c) : (c, _inclusion_score)
                  for c in self.output_labels
                  }

        data_grouped = data_tmp.groupby(
            [x for x in self.input_data.columns]).agg(
            n=(self.case_col, 'count'),
            Cases=(self.case_col, _concatenate_strings),
            **params)
        data_grouped = data_grouped[data_grouped['n'] >= self.n_cut]
        inc_columns = ['Inc_{}'.format(i) for i in self.output_labels]
        if (self.inc_score2 is None):
            data_grouped[self.output_labels] = (
                    data_grouped[inc_columns] >= self.inc_score1
            ).astype(int)
        else:
            if (self.U is None):
                raise Exception('When inc.score2 is specified,'
                                + 'U must be specified as well.')
            if (self.U != 0 and self.U != 1):
                raise Exception('U must be 0 or 1.')
            if (self.U == 1):
                data_grouped[self.output_labels] = (data_grouped[inc_columns] >
                                                    self.inc_score1).astype(int)
            if (self.U == 0):
                data_grouped[self.output_labels] = (data_grouped[inc_columns] >=
                                                    self.inc_score2).astype(int)

        res = data_grouped.reset_index()
        if self.rename_columns:
            rename_dic = {
                k: v for k, v in
                zip(self.input_labels, COLUMN_LABELS[:len(self.input_labels)])
            }
            self.rename_dictionary = rename_dic
            res.columns = map(lambda x: rename_dic[x] if x in rename_dic.keys()
            else x, res.columns)
            l = len(self.input_labels)
            self.input_labels = COLUMN_LABELS[:l]
        convert_dict = dict()
        for col in inc_columns:
            convert_dict[col]='float64'
        for col in self.output_labels:
            convert_dict[col]='int64'
        self.preprocessed_data_raw = res.astype(convert_dict)
        self.preprocessed_data = res[[x for x in self.input_labels]
                                     + self.output_labels]
        self.preprocessing = True

    def get_preprocessed_data(self, raw=False):
        '''
        Function to derive the truth table with the provided parameters
        in OptimizationContext.
        Returns
        -------
        preprocessed_data = dataframe
                  Truth table derived from the data.
        Examples
        --------
        >>> df = pd.DataFrame([[1,0,0,1],[0,1,0,1],[1,1,1,0],[1,1,1,1],[1,1,1,1]],
        ...                  columns = ["A","B","C","O"])
        >>> context = OptimizationContext(df,["O"],inc_score1=0.5)
        >>> context.get_preprocessed_data()
            A  B  C  O
        0   0  1  0  1
        1   1  0  0  1
        2   1  1  1  1
        '''
        if not self.preprocessing:
            self._preprocess_data()
        if raw:
            return self.preprocessed_data_raw
        else:
            return self.preprocessed_data

    # Function to derive the truth table from a data frame
    def _get_levels(self):
        inputs = self.preprocessed_data.drop(self.output_labels, axis=1)
        if len(self.input_labels) == 1:
            dim_corrected = [inputs.iloc[:, 0].values]
        else:
            dim = [inputs[col].unique() for col in inputs]

            dim_corrected = []
            for ar in dim:
                if len(ar) > 1:
                    dim_corrected.append(ar)
                else:
                    dim_corrected.append([0, 1])
        levels = [len(x) for x in dim_corrected]

        return levels


    def _prepareRows(self):
        if not self.preprocessing:
            self._preprocess_data()
        mask1 = self.preprocessed_data[self.output_labels]
        nr_rows = len(self.preprocessed_data.index)
        n = len(mask1.columns)
        non_zero_output = [1] * n
        multi_mask = self.preprocessed_data[
            self.output_labels].isin(non_zero_output)
        mask = multi_mask.aggregate(any, axis=1)
        positiveRows = self.preprocessed_data[mask]
        columns = [col for col in positiveRows.columns if
                   col not in self.output_labels]
        positiveInputs = positiveRows[columns]
        positiveInputs_rownames = list(positiveInputs.index)
        inputs = self.preprocessed_data.drop(self.output_labels, axis=1)

        if len(self.input_labels) == 1:
            dim_corrected = [inputs.iloc[:, 0].values]
        else:
            dim = [pd.unique(col.values).tolist() for _, col in
                   inputs.iteritems()]
            dim_corrected = []
            for ar in dim:
                if len(ar) > 1:
                    dim_corrected.append(ar)
                else:
                    dim_corrected.append([0, 1])
        levels = [len(x) for x in dim_corrected]
        cares_indexes = list()

        allInputs = pd.DataFrame(itertools.product(*dim_corrected))
        zero_output = [0] * n
        multi_mask_zero = self.preprocessed_data[
                self.output_labels].isin(zero_output)
        mask_zero = multi_mask_zero.aggregate(all, axis=1)
        negativeRows = self.preprocessed_data[mask_zero]
        negativeInputs = negativeRows[columns]
        indexes = list()
        for x in negativeInputs.values:
                ind = _find_index2(allInputs, x)
                indexes.append(ind[0])

        allInputs_table = allInputs.drop(indexes)
        allInputs_table.columns = columns

        for x in positiveInputs.values:
            ind = _find_index2(allInputs_table, x)
            cares_indexes.append(ind[0])

        if self.multi_output:
            outcols_merge = pd.merge(allInputs_table, positiveRows,
                                         how='inner', on=columns)
            outcols = outcols_merge[self.output_labels]
            self.outputcolumns = outcols.transpose().to_numpy()
        else:
            self.outputcolumns = [1] * nr_rows

        self.levels = levels
        self.cares = cares_indexes
        self.table = allInputs_table.to_numpy()
        self.labels = columns
        self.positive_cares = positiveInputs_rownames

        self.prepare_rows_called = True


    def _get_prime_implicants_ON_DC(self):

        if not self.prepare_rows_called:
            self._prepareRows()
        if len(self.table) == 0:
            prime_implicants = tuple()
            return prime_implicants

        table = self.table.astype(int).tolist()
        column_number = len(table[0])
        preprocessed_table = _preprocess_input(table)
        prime_implicants = []
        groups = _create_groups(preprocessed_table, column_number,
                                self.cares, self.outputcolumns,
                                self.multi_output)
        reduction_nr = 0
        while (True):
            reduction_res = _reduction_step(groups, column_number,
                                            self.multi_output)
            if ((reduction_res)['implicants']):
                for i in reduction_res['implicants']:
                    prime_implicants.append(i)
            if (reduction_res['reduction'] is False):
                break
            groups = reduction_res['groups']
            reduction_nr += 1

        prime_implicants = _eliminate_minterms(table,
                                               prime_implicants,
                                               self.levels,
                                               self.multi_output)
        coverage_dict = _create_care_translation_dict(self.cares,
                                                      self.positive_cares)


        if self.multi_output:

            prime_implicants_fin = tuple(ImplicantMultiOutput(
                self,
                _minterm_to_str(x[0],
                                self.levels,
                                self.labels,
                                0,
                                self.multi_output
                                ),
                x[0],
                {coverage_dict[y] for y in x[1]},
                outputs=list(x for x in x[2]),
                output_labels=[self.output_labels[i - 1]
                               for i in list(x for x in x[2])])
                                         for x in prime_implicants
                                         )
        else:
            essential_indexes = _calculate_essential_indexes(prime_implicants)
            if len(essential_indexes) > 0:
                essential_implicants = [x for x in prime_implicants if
                                        set(x[1]).intersection(
                                            essential_indexes)]

                all_essential_indexes = set.union(*[set(x[1])
                                                    for x in
                                                    essential_implicants])
                useless = [x for x in prime_implicants
                           if x[1].issubset(all_essential_indexes) and
                           x not in essential_implicants]
                prime_implicants_new = [x for x in prime_implicants
                                        if x not in useless]
                essential_tag = [x in essential_implicants
                                 for x in prime_implicants_new]
                prime_implicants_fin = tuple(Implicant(
                    self,
                    _minterm_to_str(x[0],
                                    self.levels,
                                    self.labels,
                                    0,
                                    self.multi_output),
                    x[0],
                    {coverage_dict[y] for y in x[1]}, essential=y)
                                             for x, y in
                                             zip(prime_implicants_new,
                                                 essential_tag)

                                             )
            else:
                prime_implicants_fin = tuple(Implicant(
                    self,
                    _minterm_to_str(x[0],
                                    self.levels,
                                    self.labels,
                                    0,
                                    self.multi_output),
                    x[0],
                    {coverage_dict[y] for y in x[1]})
                                             for x in prime_implicants

                                             )
        return prime_implicants_fin


    def _output_coverage_of_pi(self, raw_implicant):
        data = self.preprocessed_data
        outputs = self.output_labels
        res = set()
        for ind, out in enumerate(outputs):
            if all(data[out][data.apply(
                    lambda row_series: all(x in y for x, y in
                                           zip(row_series.values,
                                               raw_implicant)), axis=1)]):
                res.add(ind + 1)
        return res

    def _get_prime_implicants_on_off(self):

        if not self.preprocessing:
            self._preprocess_data()
        #  constant outputs
        outputs = self.preprocessed_data[self.output_labels]

        if (all(outputs.apply(lambda col: True if len(col.unique()) == 1
        else False, axis=0))):
            return self._get_prime_implicants_ON_DC()

        if len(self.output_labels) > 1:
            self.levels = self._get_levels()
            self.labels = [col for col in self.preprocessed_data.columns if
                           col not in self.output_labels]
            self.cares = [int(x) for x in
                          self.preprocessed_data[
                              self.output_labels].apply(
                              lambda row: any(x for x in row), axis=1).index]
            onset, offset = on_off_grouping_mo(self.preprocessed_data,
                                               self.output_labels)

            impl_dict = reduction_mo(onset, offset)
            tmp_res = []
            for tag, implicants in impl_dict.items():

                for im in implicants:
                    raw_im = _transform_to_raw_implicant(im, self.levels)
                    o_tag = self._output_coverage_of_pi(raw_im)

                    cov = frozenset([int(x) for x in
                                     self.preprocessed_data[
                                         self.preprocessed_data.apply(
                                             lambda row: all(x in y for x, y in
                                                             zip(row, raw_im)),
                                             axis=1)].index])

                    tmp_res.append(ImplicantMultiOutput(self,
                                                        _minterm_to_str(raw_im,
                                                                        self.levels,
                                                                        self.labels,
                                                                        0,
                                                                        self.multi_output),
                                                        raw_im,
                                                        cov,
                                                        o_tag,
                                                        self.output_labels))
            res = []
            for i in tmp_res:
                add_to_res = True
                for x in res:
                    if x.raw_implicant == i.raw_implicant:
                        x.outputs.update(i.outputs)
                        add_to_res = False
                        break
                if add_to_res:
                    res.append(i)

            for x in res:
                x.outputs = list(x.outputs)

            return res



        else:
            if all(self.preprocessed_data[self.output_labels[0]]):
                return self._get_prime_implicants_ON_DC()
            self.levels = self._get_levels()

            self.labels = [col for col in self.preprocessed_data.columns if
                           col not in self.output_labels]

            self.cares = [int(x) for x in self.preprocessed_data[
                self.preprocessed_data[self.output_labels[0]] == 1].index]

            prime_implicants = []
            essentials = get_essential_implicants(self.preprocessed_data,
                                                  self.output_labels[0],
                                                  self.levels)
            if len(essentials) > 0:
                cov_essentials = _reduce_the_onset(essentials,
                                                   self.preprocessed_data,
                                                   self.output_labels[0])
                for impl, cov in zip(essentials, cov_essentials):
                    raw_i = _transform_to_raw_impl(impl, self.levels)

                    prime_implicants.append(Implicant(
                        self,
                        _minterm_to_str(raw_i,
                                        self.levels,
                                        self.labels,
                                        0,
                                        self.multi_output),
                        raw_i, cov, essential=True))
                    reduced_indexes = set.union(
                        *[set(x) for x in cov_essentials]) if len(
                        cov_essentials) > 0 else set()
                    data_reduced = self.preprocessed_data.drop(reduced_indexes,
                                                               inplace=False)
                    onset, offset = _on_off_grouping(data_reduced,
                                                     self.output_labels[0],
                                                     self.multi_output)
            else:
                onset, offset = _on_off_grouping(self.preprocessed_data,
                                                 self.output_labels[0],
                                                 self.multi_output)
            impl_dict = _reduction(onset, offset)

            for impl, cov in impl_dict.items():
                raw_i = _transform_to_raw_implicant(impl, self.levels)

                prime_implicants.append(Implicant(
                    self,
                    _minterm_to_str(raw_i,
                                    self.levels,
                                    self.labels,
                                    0,
                                    self.multi_output),
                    raw_i, cov))

            return prime_implicants
    def _get_prime_implicants_boom(self):
        if not self.preprocessing:
            self._preprocess_data()
        pass
        #while(100):
         #  candidates =  best_literals(self.preprocessed_data,
          #                            self.output_labels,
           #                           self.case_col)



    def get_prime_implicants(self):
        """
        Function computes the prime implicants.
        Returns
        -------
        prime_implicants = tuple
        A tuple of prime implicant objects.
        Notes
        -------
        The string representation of a prime implicant object is such that:
        I.) Binary case:
        a) a positive literal is printed in upper case,
        b) a negative literal is printed in lower case, and
        c) an essential prime implicant is prefixed with a hash tag (#).
        II.) Multi-value case:
        a) all literals are coded with upper case letters
        b) the value of the literal is written in curly ({}) brackets.
        c) an essential prime implicant is prefixed with a hash tag (#).
        Examples
        --------
        >>> df = pd.DataFrame([[1,1,0,1],[0,0,1,1],[1,0,1,0],[0,1,0,1]],
        ...                   columns=["A","B","C","OUT"])
        >>> context = cora.OptimizationContext(data = df,
        ...                                    output_labels = ["OUT"])
        >>> context.get_prime_implicants()
        (B, c, #a)
        >>> df = pd.DataFrame([[1,2,0,1,1,2,1],
        ...                   [1,1,1,0,2,0,0],
        ...                   [0,2,1,0,0,1,2],
        ...                   [0,2,2,0,1,1,1],],
        ...                  columns=["A","B","C","D","OUT1","OUT2","OUT3"])
        >>> context = cora.OptimizationContext(data = df,
        ...                                    output_labels =["OUT1{1,2}",
        ...                                                    "OUT2{1}",
        ...                                                     "OUT3{1,0}"],
        ...                                     algorithm = 'ON-OFF')
        >>> context.get_prime_implicants()
        [C{2}, A{1}, B{1}, D{1}, C{0}, B{2}*C{1}, A{0}, B{2}*D{0}]
        """

        if self.prime_implicants is not None:
            return self.prime_implicants

        if self.algorithm == "Boom":
            self.prime_implicants = self._get_prime_implicants_boom()
        elif self.algorithm == "ON-DC":
            self.prime_implicants = self._get_prime_implicants_ON_DC()
        elif self.algorithm == "ON-OFF":
            self.prime_implicants = self._get_prime_implicants_on_off()
        else:
            raise AttributeError(
                'Unknown algorithm "{}"?'.format(self.algorithm))

        return self.prime_implicants


    def prime_implicant_chart(self):
        """
        Function creates a prime implicant chart (coverage matrix).
        Returns
        -------
        dataframe
        columns: The labels of the columns correspond to the truth table rows
                 with positive outputs.
        rows: prime implicants
        If a prime implicant covers a positive row of the truth table, it
        receives a 1 entry in the corresponding position; an entry 0 otherwise.
        Example
        -------
        >>> df = pd.DataFrame([[1,1,0,1],
        ...           [0,0,1,1],
        ...           [1,0,1,0],
        ...           [0,1,0,1]],
        ...           columns=["A","B","C","OUT"])
        >>> context = cora.OptimizationContext(data = df,
        ...                                    output_labels = ["OUT"])
        >>> context.prime_implicant_chart()
            0  1  3
        B   0  1  1
        c   0  1  1
        #a  1  1  0
        """
        if self.pi_chart is not None:
            return self.pi_chart
        prime_implicants = self.get_prime_implicants()
        cares = set().union(*[set(x.coverage) for x
                              in prime_implicants])

        res = np.zeros((len(cares), len(prime_implicants)), dtype=bool)
        res_idx_to_care = {v: k for k, v in enumerate(cares)}
        for row_nr, implicant in enumerate(prime_implicants):
            for x in implicant.coverage:
                if x not in cares:
                    continue
                column_nr = res_idx_to_care[x]
                res[column_nr, row_nr] = True
        if self.multi_output:
            return pd.DataFrame(res.transpose(),
                                columns = cares,
                                index = ['{}, {}'.format(x.implicant, x.outputs)
                                       for x in prime_implicants]).astype(int)
        self.pi_chart = pd.DataFrame(
            res.transpose(),
            columns = cares,
            index = [(x.implicant) for x in prime_implicants]).astype(int)
        return self.pi_chart

    def get_irredundant_sums(self, max_depth=None):
        """
        Function computes all irredundat sums.
        Parameters
        ----------
        max_depth : int
                   A positive integer denoting max number of prime implicants
                   in the solution.
        Returns
        -------
        irredundat_objects : array of objects
        The final array containing all irredundant sums
        with respect to the max_depth condition.
        Notes
        ------
        Note that the number of all irredundant sums can grow exponentially.
        Requesting all irredundant sums for a large prime implicant chart
        might consume a lot of memory and lead to out-of-memory crashes.
        Examples
        --------
        >>> df = pd.DataFrame([[1,1,0,1],
        ...           [0,0,1,1],
        ...           [1,0,1,0],
        ...           [0,1,0,1]],
        ...           columns=["A","B","C","OUT"])
        >>> context = cora.OptimizationContext(data = df,
        ...                                    output_labels = ["OUT"])
        >>> context.get_irredundant_sums()
        [M1: #a + B, M2: #a + c]
        """
        if self.irredundat_sums is not None:
            return self.irredundat_sums
        prime_implicants = self.get_prime_implicants()
        if len(prime_implicants) == 0:
            return []

        if self.multi_output:
            raise RuntimeError("irredudant sums are not supported in multi\
                         output mode. Use get_irredundant_systems")

        result = _find_irredundant_sums_native(([(i, i.coverage)
                                                 for i in prime_implicants]),
                                               self.cares

                                               )
        irredundant_objects = []
        for i, system in enumerate(result):
            irredundant_objects.append(IrredundantSystem(
                self, system, i + 1))

        self.irredundat_sums = irredundant_objects
        return self.irredundat_sums




    def system_details(self):
        """
        Function gives the statistical overview of a solution system.
        Returns
        -------
        df_final : dataframe
        (a) system unique identification label
        (b) systems' string representation (for one-output systems only)
        (c) the coverage and the inclusion score of the system
        Examples
        --------
        >>> df = pd.DataFrame([[1,1,0,1],
        ...           [0,0,1,1],
        ...           [1,0,1,0],
        ...           [0,1,0,1]],
        ...           columns=["A","B","C","OUT"])
        >>> context = cora.OptimizationContext(data = df,
        ...                                    output_labels = ["OUT"])
        >>> context.system_detials()
                          Cov.  Inc.
        Solution details   1.0   1.0
        """
        if self.sol_details is not None:
            return self.sol_details
        if not self.multi_output:
            irr_sums = self.get_irredundant_sums()
            solution_sample = irr_sums[0]

        else:
            irr_systems = self.get_irredundant_systems()
            solution_sample = irr_systems[0]

        self.sol_details = pd.DataFrame({
            'Cov.': round(solution_sample.coverage_score(), 2),
            'Inc.': round(solution_sample.inclusion_score(), 2)},
            index=["Solution details"])

        return self.sol_details



    def pi_details(self):
        """
        Function gives the statistical overview of all prime implicants.
        Returns
        -------
        df_implicant : dataframe
        Dataframe collects all statistical values of an implicant.
        Note
        ------
       a) The inclusion score of a prime implicant is the ratio between the
       number of cases that are covered by this prime implicant and show the
       analyzed outcome, and the number of cases that are covered by this
       prime implicant.
       b) The coverage score of a prime implicant is the ratio between the
       number of cases that are covered by this prime implicant and show the
       analyzed outcome, and the number of cases that show the analyzed outcome.
        Examples
        --------
        >>> df = pd.DataFrame([[1,1,0,1],
        ...           [0,0,1,1],
        ...           [1,0,1,0],
        ...           [0,1,0,1]],
        ...           columns=["A","B","C","OUT"])
        >>> context = cora.OptimizationContext(data = df,
        ...                                    output_labels = ["OUT"])
        >>> context.pi_detials()
           PI  Cov.r  Inc.    M1    M2
        0   B   0.67   1.0  0.33   NaN
        1   c   0.67   1.0   NaN  0.33
        2  #a   0.67   1.0  0.33  0.33
        """
        if self.details is not None:
            return self.details
        prime_implicants = self.get_prime_implicants()

        cov_x = [(x.implicant,
                  round(x.coverage_score(), 2),
                  round(x.inclusion_score(), 2)) for x in prime_implicants]

        if not self.multi_output:
            irr_sums = self.get_irredundant_sums()
            new_cols = ["M" + str(x.index) for x in irr_sums]
            cov_per_impl = [x.impl_cov_score() for x in irr_sums]

        else:
            irr_systems = self.get_irredundant_systems()
            new_cols = ["S" + str(x.index) for x in irr_systems]
            cov_per_impl = [x.impl_cov_score() for x in irr_systems]

        df_implicant = pd.DataFrame(cov_x, columns=['PI', 'Cov.r', 'Inc.'])
        for f_cov_per_impl, f_label in zip(cov_per_impl, new_cols):
            tmp = [round(f_cov_per_impl[x.implicant], 2) if x.implicant
                                            in f_cov_per_impl.keys() else None
                   for x in prime_implicants]
            df_implicant[f_label] = tmp

        self.details = df_implicant
        return self.details


    def get_irredundant_systems(self):
        """
        Function computes all irredundat systems for input data with
        multiple output columns.
        Returns
        -------
        res : array of objects
        (a) unique identification label of the irredundant system object,
        (b) solution systems containing solution string per output.
        Note
        ------
        Individual functions inside one system may not be irredundant,
        but the entire system is always irredundant.
        Example
        --------
        >>> df = pd.DataFrame([[1,2,0,1,1,2,1],
        ...                    [1,1,1,0,2,0,0],
        ...                    [0,2,1,0,0,1,2],
        ...                    [0,2,2,0,1,1,1],],
        ...                    columns=["A","B","C","D","OUT1","OUT2","OUT3"])
        >>> context = cora.OptimizationContext(df, ["OUT1{1,2}",
        ...                                         "OUT2{1}",
        ...                                         "OUT3{1,0}"],
        ...                                          algorithm = 'ON-OFF')
        >>> context.get_irredundant_systems()
        [---- System 1 ----
        OUT1{1, 2}: B{1}
        OUT2{1}: B{2}*D{0}
        OUT3{0, 1}: B{1}
        , ---- System 2 ----
        OUT1{1, 2}: A{1}
        OUT2{1}: B{2}*D{0}
        OUT3{0, 1}: A{1}
        ]
        """
        if self.irredundant_systems is not None:
            return self.irredundant_systems

        prime_implicants = self.get_prime_implicants()

        if not self.multi_output:
            raise RuntimeError("irredudant systems are not supported in single\
                           output mode. Use get_irredundant_sums")
        res, l = self._single_ir_systems_for_multi_output()

        mult_input = [set(frozenset(imp for imp in irs) for irs in f) for f in
                      res]
        reduction_result = reduce(_boolean_multiply, mult_input)
        res = []
        index = 0
        for r in reduction_result:
            index = index + 1
            single_res = []
            for j in range(l):
                single_res.append([prime_implicants[i] for i in r if j + 1 in
                                   prime_implicants[i].outputs])
            res.append(IrredundantSystemsMulti(self, single_res,
                                               index
                                               ))
        self.irredundant_systems = res
        return self.irredundant_systems

    def _single_ir_systems_for_multi_output(self):
        prime_implicants = self.get_prime_implicants()
        l = len(self.output_labels)
        imp_per_output = []
        for i in range(l):
            imp_per_output.append(list())
        for i, impl in enumerate(prime_implicants):
            for j in range(l):
                if (j + 1) in set(impl.outputs):
                    imp_per_output[j].append((i, impl))
        res = []
        for k, system in enumerate(imp_per_output):
            coverage = set().union(*[set(i[1].coverage) for i in system])
            result = _find_irredundant_sums_native(([(i, impl.coverage)
                                                     for i, impl in system]),
                                                   coverage)
            irredundant_objects = []
            for i, system in enumerate(result):
                irredundant_objects.append(IrredundantSystem(self,
                                                             system,
                                                             i + 1
                                                             ))

            for i in irredundant_objects:
                res.append([i.system for i in irredundant_objects])

        return res, l

    def get_solution_dataframe(self):
        """
        Return
        ------
        dataframe :
        The summary dataframe separately specifies an occurence
        of a prime implicant inside each output within a system.
        Example
        --------
        >>> df = pd.DataFrame([[1,2,0,1,1,2,1],
        ...                    [1,1,1,0,2,0,0],
        ...                    [0,2,1,0,0,1,2],
        ...                    [0,2,2,0,1,1,1],],
        ...                    columns=["A","B","C","D","OUT1","OUT2","OUT3"])
        >>> context = cora.OptimizationContext(df, ["OUT1{1,2}",
        ...                                         "OUT2{1}",
        ...                                         "OUT3{1,0}"],
        ...                                          algorithm = 'ON-OFF')
        >>> context.get_solution_dataframe()
           C{2}  A{1}  B{1}  D{1}  C{0}  B{2}*C{1}  A{0}  B{2}*D{0} Output System
        0     0     0     1     0     0          0     0          0      1      1
        1     0     0     0     0     0          0     0          1      2      1
        2     0     0     1     0     0          0     0          0      3      1
        3     0     1     0     0     0          0     0          0      1      2
        4     0     0     0     0     0          0     0          1      2      2
        5     0     1     0     0     0          0     0          0      3      2
        """
        prime_implicants = self.get_prime_implicants()
        if not self.multi_output:

            solutions = self.get_irredundant_sums()

            n_rows = len(solutions)
            n_cols = len(prime_implicants)

            data = np.zeros(shape=(n_rows, n_cols), dtype=int)
            for r, solution in enumerate(solutions):
                impl_set = frozenset(i.implicant for i in solution.system)
                for c, implicant in enumerate(prime_implicants):
                    if implicant.implicant in impl_set:
                        data[r, c] = 1
            self.solution_dataframe = pd.DataFrame(data, range(n_rows),
                                                   [i.implicant for i in
                                                    prime_implicants])
            return self.solution_dataframe

        else:

            solutions = self.get_irredundant_systems()
            l = len(self.output_labels)
            n_rows = len(solutions) * l
            n_cols = len(prime_implicants)
            data = np.zeros(shape=(n_rows, n_cols), dtype=int)
            for i, solution in enumerate(solutions):

                for out, single_output_solution in enumerate(
                        solution.system_multiple):
                    impl_set_out = frozenset(
                        i.implicant for i in single_output_solution)
                    for c, implicant in enumerate(prime_implicants):
                        if implicant.implicant in impl_set_out:
                            data[l * i + out, c] = 1
            df = pd.DataFrame(data, range(n_rows),
                              [i.implicant for i in prime_implicants])

            #out_strings = [str(i + 1) for i in range(l)]
            out_strings=self.output_labels_final
            df["Output"] = out_strings * len(solutions)
            tmp = []
            for i in range(len(solutions)):
                tmp.extend([str(i + 1)] * l)
            df["System"] = tmp
            return df

    def _get_network_representation(self):
        if self.solution_dataframe is None:
            self.get_solution_dataframe()

        n = self.solution_dataframe.shape[0]
        res = pd.DataFrame(np.zeros((50, 50)), dtype=int)
        for i in range(0, 50):
            for j in range(i + 1, 50):
                x = self.solution_dataframe.iloc[[i, j]].astype(bool).apply(
                    lambda x: x[i] ^ x[j], axis=0).astype(int).sum()
                res.at[i, j] = x
                res.at[j, i] = x
        return res


"""
 Class represents a single irredundant system.
 
The complete solution space of multi-output minimization is the set of 
all systems.
 
Parameters
----------
system : array of strings
         An array of prime implicants forming the solution.  
index : int
        An index of the solution.
cov_score : float
            statistical value
incl_score : float
             statistical value
 
"""
class IrredundantSystemsMulti():
    # Class where self represents just one single system
    # the whole solution is composed of these systems

    def __init__(self, context, system_multiple, index
                 ):
        self.context = context
        self.system_multiple = system_multiple
        self.index = index
        self.output_labels = context.output_labels_final
        self.cov_score = None
        self.incl_score = None

    def __str__(self):

        res = ""
        res += '---- System {} ----\n'.format(self.index)
        for j, system in enumerate(self.system_multiple):
            if any(str(impl.implicant) == '1' for impl in system):
                res += ('{}: 1\n'.format(self.output_labels[j]))
            elif system == []:
                res += ('{}: 0\n'.format(self.output_labels[j]))
            else:
                res += ('{}: {}\n'.format(self.output_labels[j],
                                          ' + '.join(impl.implicant
                                                     for impl in system)))
        res += '\n'
        return res

    def get_descriptive_string(self, cov):
        res = ""
        res += '---- System {} ----\n'.format(self.index)
        for j, system in enumerate(self.system_multiple):

            if any(str(impl.implicant) == '1' for impl in system):
                res += ('1 <=> {}\n'.format(self.output_labels[j]))
            elif system == []:
                res += ('0 <=> {}\n'.format(self.output_labels[j]))
            else:
                if (self.context.inc_score2 is not None and
                        self.context.U == 0):
                    final_inc_score = self.context.inc_score2
                else:
                    final_inc_score = self.context.inc_score1

                solution_cov = self.coverage_score()
                solution_inc = self.inclusion_score()
                if (solution_inc >= final_inc_score
                        and solution_inc >= 0.5
                        and solution_cov >= cov
                        and solution_cov >= 0.5):
                    res += ('{1} <=> {0}\n'.format(self.output_labels[j],
                                                   ' + '.join(impl.implicant
                                                              for impl in
                                                              system)))

                elif (solution_inc >= final_inc_score
                      and solution_inc >= 0.50):
                    res += ('{1} => {0}\n'.format(self.output_labels[j],
                                                  ' + '.join(impl.implicant
                                                             for impl in
                                                             system)))
                elif (solution_cov >= cov and solution_cov >= 0.50):
                    res += ('{1} <= {0}\n'.format(self.output_labels[j],
                                                  ' + '.join(impl.implicant
                                                             for impl in
                                                             system)))
                else:

                    res = "Warning! "
                    break
        return res

    def __repr__(self):
        return str(self)

    def unique_implicants(self):
        return [x for x in set(sum(self.system_multiple, list()))]

    def sort_implicants(self):
        implicants = self.unique_implicants()
        res = defaultdict(list)

        for imp1 in implicants:
            for imp2 in implicants:
                if (len(imp1.outputs) == 1
                        and set(imp1.outputs).issubset(imp2.outputs)
                        and imp1 != imp2):
                    res[imp1].append(imp2)

                elif (len(imp1.outputs) > 1
                      and set(imp1.outputs) == set(imp2.outputs)
                      and imp1 != imp2):
                    res[imp1].append(imp2)
        return res

    def impl_cov_score(self):
        
        '''
        Return:
        Float value, representing the proportion of instantiations of the implicant among
        all instantiations of the (corresponding) outcomes, combinations of outcomes, respectively, 
        unique to the implicant.
        '''
            
        data = self.context.data
        input_columns = self.context.input_data.columns
        output_columns = self.context.output_labels
        implicants = self.unique_implicants()
        sorted_implicants = self.sort_implicants()

        unique_cov = []

        for impl in sorted_implicants:
            tmp_data = data[input_columns]
            outputs = dict()
            for indx,output in enumerate(output_columns):
                if int(indx+1) in impl.outputs:
                    outputs[output] = 1
                else:
                    outputs[output] = 0
            mask = data.apply(lambda row_series: all(row_series[key] == value
                                                     for key , value
                                                     in outputs.items()),

                              axis=1)

            tmp_coresponding_data = tmp_data[mask]
            tmp = tmp_coresponding_data.apply(
                lambda row_series: row_series.name if all(x in y for x, y in
                                                          zip(row_series.values,
                                                            impl.raw_implicant))
                else np.NAN, axis=1)

            s_in = set(x for x in tmp.values if not np.isnan(x))
            s_out = set()
            if len(sorted_implicants[impl]) > 0:
                for imp_2 in sorted_implicants[impl]:
                    tmp2 = tmp_coresponding_data.apply(
                        lambda row_series: row_series.name if all(
                            x in y for x, y in
                            zip(row_series.values,
                                imp_2.raw_implicant))
                        else np.NAN, axis=1)
                    s_out.update(set(x for x in tmp2.values if not np.isnan(x)))


            unique_cov.append(len(s_in - s_out) / len(tmp_coresponding_data.index))
        return {str(impl_i.implicant): coverage
                for impl_i, coverage in zip(sorted_implicants, unique_cov)}

        # coverage of the system

    def coverage_score(self):
        
        '''
        Return:
        Float value, representing the proportion of instantiations of the implicant among
        all instantiations of the (corresponding) outcomes, combinations of outcomes, respectively.
        '''
        
        if self.cov_score is None:
            data = self.context.data
            input_columns = self.context.input_data.columns
            output_columns = self.context.output_labels
            implicants = self.unique_implicants()

            tmp_data = data[data.apply(lambda row_series:
                                       any(row_series[output] == 1
                                           for output in output_columns),
                                       axis=1)][
                input_columns]

            self.cov_score = tmp_data.apply(
                lambda row_series: 1.0 if any(all(x in y for x, y in
                                                  zip(row_series.values,
                                                      i.raw_implicant))
                                              for i in implicants)
                else 0.0, axis=1).mean()
        return self.cov_score

    # inclusion of the system
    def inclusion_score(self):
        '''
        Return:
        Float value, representing the proportion of instantiations of the disjunction of implicants
        in conjunction with the outcomes among all instantiations of the disjunction of implicants.
        '''
        if self.incl_score is None:
            data = self.context.data
            input_columns = self.context.input_data.columns
            output_columns = self.context.output_labels
            tmp_data = data[input_columns]
            implicants = self.unique_implicants()

            mask = tmp_data.apply(
                lambda row_series: any(all(x in y for x, y in
                                           zip(row_series[input_columns].values,
                                               i.raw_implicant))
                                       for i in implicants), axis=1)
            new_out = data.apply(
                lambda row_series: (
                    1 if sum(row_series[output_columns].values) > 0
                    else 0), axis=1)

            data["new_output"] = new_out
            self.incl_score = data.loc[mask, "new_output"].mean()
        return self.incl_score

class IrredundantSystem():
    """
    Class represents an irredundant system of an optimization problem with a single outcome.
     
    Parameters
    ----------
    system : array of strings
             An array of implicants which form the solution.
    index : int
            An index of the solution.
    cov_score : float
                statistical value
    incl_score : float
                 statistical value

    """

    def __init__(self, context, system, index):
        self.context = context
        if len(context.output_labels) == 1:
            self.system = sorted(system, key=lambda x: x.essential,
                                 reverse=True)
        else:
            self.system = system
        self.index = index
        self.cov_score = None
        self.incl_score = None
        self.output = context.output_labels_final[0]

    def __str__(self):
        return 'M{}: {}'.format(self.index,
                                ' + '.join(str(i.implicant)
                                           for i in self.system))

    def get_descriptive_string(self, cov):
        solution_cov = self.coverage_score()
        solution_inc = self.inclusion_score()

        if (self.context.inc_score2 is not None and
                self.context.U == 0):
            final_inc_score = self.context.inc_score2
        else:
            final_inc_score = self.context.inc_score1

        if (solution_cov >= cov
                and solution_cov >= 0.5
                and solution_inc >= 0.5
                and solution_inc >= final_inc_score):
            return '{} <=> {}'.format(' + '.join(str(i.implicant)
                                                 for i in self.system),
                                      self.output)
        elif (solution_cov >= cov and solution_cov >= 0.5):
            return '{} <= {}'.format(' + '.join(str(i.implicant)
                                                for i in self.system),
                                     self.output)
        elif (solution_inc >= final_inc_score and solution_inc >= 0.5):
            return '{} => {}'.format(' + '.join(str(i.implicant)
                                                for i in self.system),
                                     self.output)
        else:
            return "Warning!"

    def impl_cov_score(self):
        
        '''
        Return:
        Float value, representing the proportion of instantiations of the implicant among
        all instantiations of the outcome unique to the implicant.
        '''
            
        data = self.context.data
        input_columns = self.context.input_data.columns
        output_column = self.context.output_labels[0]

        tmp_data = data[input_columns]
        data[output_column] == 1
        tmp_positive_data = tmp_data[data[output_column] == 1]
        impl_cov = []
        cov_count = {}
        for i, impl_i in enumerate(self.system):

            tmp = tmp_positive_data.apply(
                lambda row_series: row_series.name if all(x in y for x, y in
                                                        zip(row_series.values,
                                                        impl_i.raw_implicant))
                else None, axis=1)
            s = set(x for x in tmp.values if not np.isnan(x))
            impl_cov.append(s)
            for x in s:
                if x in cov_count.keys():
                    cov_count[x] += 1
                else:
                    cov_count[x] = 1
        unique_cov = [{x for x in ic if cov_count[x] == 1} for ic in impl_cov]

        return {str(impl_i.implicant): len(ic) / len(tmp_positive_data.index)
                for impl_i, ic in zip(self.system, unique_cov)}

    def __repr__(self):
        return str(self)

    def nr_implicants(self):
        return len(self.system)

    def coverage_score(self):
        
        '''
        Return:
        Float value, representing the proportion of instantiations of the implicant among
        all instantiations of the outcome.
        '''
        
        data = self.context.data
        input_columns = self.context.input_data.columns
        output_column = self.context.output_labels[0]

        tmp_data = data[data[output_column] == 1][input_columns]
        self.cov_score = tmp_data.apply(
            lambda row_series: 1.0 if any(all(x in y for x, y in
                                              zip(row_series.values,
                                                  i.raw_implicant))
                                          for i in self.system)
            else 0.0, axis=1).mean()
        return self.cov_score

    def inclusion_score(self):
    
        '''
        Return:
        Float value, representing the proportion of instantiations of the disjunction of implicants
        in conjunction with the outcome among all instantiations of the disjunction of implicants.
        '''
        
        data = self.context.data
        input_columns = self.context.input_data.columns
        output_column = self.context.output_labels[0]
        tmp_data = data[input_columns]
        mask = tmp_data.apply(
            lambda row_series: any(all(x in y for x, y in
                                       zip(row_series[input_columns].values,
                                           i.raw_implicant))
                                   for i in self.system), axis=1)
        self.incl_score = data.loc[mask, output_column].mean()
        return self.incl_score


# Class which defines a multi-output minterm/item and describes its properties
# (coverage, is_reduced, tag) over the reduction process.
#
# Parameters
# ----------
#
#    minterm : A tuple of numbers. In the first iteration, every tuple
#              correponds to one row from the data.

#    coverage : Set of the indexes of all the rows which are covered by the
#               minterm. At the beggining of the  minimalization coverage, it
#               contains just a single number - the index of a row represented
#               by the minterm.

#    is_reduced : True when at least one reduction was performed.
#
#    tag : An arbitrary number, representing the outputs.



class MultipleOutputMinterm:

    def __init__(self, minterm, coverage, tag):
        self.minterm = tuple(x for x in minterm)
        self.coverage = frozenset(x for x in coverage)
        self.is_reduced = False
        self.tag = frozenset(x for x in tag)

    def __str__(self):
        return ('{0}, {1},{2}, tag={3}'.format(
            self.minterm,
            self.coverage,
            self.is_reduced,
            self.tag))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (self.minterm == other.minterm and
                self.coverage == other.coverage and
                self.is_reduced == other.is_reduced and
                self.tag == other.tag)

    def __hash__(self):
        return hash((self.minterm, self.coverage, self.is_reduced))

    def can_be_reduced(self, other):
        if len(self.tag.intersection(other.tag)) == 0:
            return False
        n = len(self.minterm)
        diff = 0
        for i in range(0, n):
            if self.minterm[i] != other.minterm[i]:
                diff = diff + 1
        if diff == 1:
            return (True)
        else:
            return (False)

    def reduce(self, other):
        if not self.can_be_reduced(other):
            return None
        new_tag = self.tag.intersection(other.tag)
        if new_tag == self.tag:
            self.is_reduced = True
        if new_tag == other.tag:
            other.is_reduced = True
        n = len(self.minterm)
        new_minterm = list()
        for i in range(0, n):
            new_minterm.append({})
        for i in range(0, n):
            if (self.minterm[i] == other.minterm[i]):
                new_minterm[i] = self.minterm[i]
            else:
                new_minterm[i] = self.minterm[i].union(other.minterm[i])
        return MultipleOutputMinterm(new_minterm,
                                     self.coverage.union(other.coverage),
                                     new_tag)


# Class which defines a multi-value minterm/item and describes its properties
# (coverage, is_reduced) over the reduction process.
#
# Parameters
# ----------
#
#    minterm : A tuple of numbers. In the first itteration, every tuple
#              correponds to one row from the data.

#    coverage : Set of the indexes of all the rows which are covered by the
#               minterm. At the beginning of the minimalization the coverage
#               contains just a single number - the index of a row represented
#               by the minterm.

#    is_reduced : True when at least one reduction was performed.


class MultiValueMinterm:

    def __init__(self, minterm, coverage):
        self.minterm = tuple(x for x in minterm)
        self.coverage = frozenset(x for x in coverage)
        self.is_reduced = False

    def __str__(self):
        return ('{0}, {1},{2}'.format(
            self.minterm,
            self.coverage,
            self.is_reduced
        ))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (self.minterm == other.minterm and
                self.coverage == other.coverage and
                self.is_reduced == other.is_reduced)

    def __hash__(self):
        return hash((self.minterm, self.coverage, self.is_reduced))

    def can_be_reduced(self, other):
        n = len(self.minterm)
        diff = 0
        for i in range(0, n):
            if self.minterm[i] != other.minterm[i]:
                diff = diff + 1
        if diff == 1:
            return (True)
        else:
            return (False)

    def reduce(self, other):
        if not self.can_be_reduced(other):
            return None
        self.is_reduced = True
        other.is_reduced = True
        n = len(self.minterm)
        new_minterm = list()
        for i in range(0, n):
            new_minterm.append({})
        for i in range(0, n):
            if (self.minterm[i] == other.minterm[i]):
                new_minterm[i] = self.minterm[i]
            else:
                new_minterm[i] = self.minterm[i].union(other.minterm[i])
        return MultiValueMinterm(new_minterm,
                                 self.coverage.union(other.coverage))


class ImplicantMultiOutput:
    """
    Class which refers to the prime implicant generated from
    multi-output data and its properties and relations to the outputs.
    Parameters
     ----------
     context : object
               An object of the OptimizationContext class
               to which the implacant refers.
               It contains the original data, original output_labels
               input_labels etc...

     implicant : array of sets of numbers
                Each set corresponds to one input variable.

     raw_implicant : an array of sets of numbers
                      An arrbitrary representation of an implicant containing
                      some additional information.

     coverage : set of numbers
                The set  of the indexes of the truth table,
                which are covered by the implicant. (Columns of the PI chart)

     outputs : array of integers
               The numbers represent the indexes of the corresponding outputs
               to the implicant.

     output_labels : array of strings
                     The array contains the output labels corresponding to the
                     implicant.

     cov_u : float
          statistical value

     incl_score : float
          statistical value

    """

    def __init__(self,
                 context,
                 implicant,
                 raw_implicant,
                 coverage,
                 outputs,
                 output_labels):
        self.context = context
        self.implicant = implicant
        self.raw_implicant = raw_implicant
        self.coverage = coverage

        self.outputs = outputs
        self.output_labels = output_labels

    def __str__(self):
        return ('{}'.format(self.implicant))

    def __repr__(self):
        return str(self)

    def coverage_score(self):
        
        '''
        Return:
        Float value, representing the proportion of instantiations of the implicant among
        all instantiations of the outcomes.
        '''
        
        data = self.context.data
        input_columns = self.context.input_data.columns

        if (len(input_columns) != len(self.raw_implicant)):
            raise RuntimeError(
                'Size of input columns ({}) does not match implicant\
                    size({})'.format(len(input_columns),
                                     len(self.raw_implicant)))

        tmp_data = data[input_columns]
        outputs_complex = dict()
        for indx,output in enumerate(self.context.output_labels):
               if int(indx+1) in self.outputs:
                    outputs_complex[output] = 1
               else:
                    outputs_complex[output] = 0
        mask = data.apply(lambda row_series: all(row_series[key] == value
                                                     for key , value
                                                    in outputs_complex.items()),

                                          axis=1)

        tmp_data = tmp_data[mask]

        self.cov_score = tmp_data.apply(
            lambda row_series: 1.0 if all(x in y for x, y in
                                              zip(row_series.values,
                                                  self.raw_implicant))
                                    else 0.0, axis=1).mean()
        return self.cov_score

    def inclusion_score(self):
        '''
        Return:
        Float value, representing the proportion of instantiations of the implicant
        in conjunction with the outcomes among all instantiations of the implicant.
        '''
        data = self.context.data
        input_columns = self.context.input_data.columns

        if (len(input_columns) != len(self.raw_implicant)):
            raise RuntimeError(
                'Size of input columns ({}) does not match implicant\
                    size({})'.format(len(input_columns),
                                     len(self.raw_implicant)))
        tmp_data = data[input_columns]
        if len(self.outputs) == 1:
            tmp_positive_data = tmp_data[data[self.output_labels[0]] == 1]

            self.incl_score = tmp_positive_data.apply(
                lambda row_series:
                1.0 if all(x in y for x, y in zip(row_series.values,
                                                  self.raw_implicant))
                else 0.0, axis=1).sum() / tmp_data.apply(
                lambda row_series:
                1.0 if all(x in y for x, y in zip(row_series.values,
                                                  self.raw_implicant))
                else 0.0, axis=1).sum()
        else:

            tmp_positive_data = tmp_data[data.apply(lambda row_series:
                                                    all(row_series[output] == 1
                                                        for output in
                                                        self.output_labels),
                                                    axis=1)]

            self.incl_score = tmp_positive_data.apply(
                lambda row_series: 1.0 if all(x in y for x, y in
                                              zip(row_series.values,
                                                  self.raw_implicant))
                else 0.0, axis=1).sum() / tmp_data.apply(
                lambda row_series: 1.0 if all(x in y for x, y in
                                              zip(row_series.values,
                                                  self.raw_implicant))

                else 0.0, axis=1).sum()
        return self.incl_score

class Implicant:
    """
     Class Implicant  refers to a prime implicant and its properties.
     Parameters
     ----------
     context : object
               An object of the  OptimizationContext class
               to which the implacant refers.
               It contains the original data, output_labels,
               input_labels etc...

     implicant : array of sets of numbers
                Each set corresponds to one input variabels.

      raw_implicant  : an array of sets of number
                       An arrbitrary representation of an implicant containing
                       some additional information.

      coverage : set of numbers
                 The set  of the indexes of the truth table,
                 which are covered by the implicant. (Columns of PI chart.)

      cov_u : float
          statistical value

      incl_score : float
          statistical value

    """

    def __init__(self,
                 context,
                 implicant,
                 raw_implicant,
                 coverage,
                 essential = False):
        self.context = context
        self.implicant = "#" + str(implicant) if essential else implicant
        self.raw_implicant = raw_implicant
        self.coverage = coverage
        self.essential = essential
        self.useless = False

    def __str__(self):
        return ('{0}'.format(self.implicant))

    def __repr__(self):
        return str(self)

    def coverage_score(self):
        '''
        Return:
        Float value, representing the proportion of instantiations of the implicant among
        all instantiations of the outcome.
        '''
        data = self.context.data
        input_columns = self.context.input_data.columns
        output_column = self.context.output_labels[0]

        if (len(input_columns) != len(self.raw_implicant)):
            raise RuntimeError(
                'Size of input columns ({}) does not match implicant\
                    size({})'.format(len(input_columns),
                                     len(self.raw_implicant)))
        tmp_data = data[data[output_column] == 1][input_columns]
        self.cov_score = tmp_data.apply(
                            lambda row_series: 1.0 if all(x in y for x, y in
                                zip(row_series.values,self.raw_implicant))
                                                else 0.0, axis=1).mean()
        return self.cov_score

    def inclusion_score(self):
        '''
        Return:
        Float value, representing the proportion of instantiations of the implicant
        in conjunction with the outcome among all instantiations of the implicant.
        '''
        data = self.context.data
        input_columns = self.context.input_data.columns
        output_column = self.context.output_labels[0]
        if (len(input_columns) != len(self.raw_implicant)):
            raise RuntimeError(
                'Size of input columns ({}) does not match implicant\
                    size({})'.format(len(input_columns),
                                     len(self.raw_implicant)))
        tmp_data = data[input_columns]
        tmp_positive_data = tmp_data[data[output_column] == 1]

        self.incl_score = tmp_positive_data.apply(
            lambda row_series:
            1.0 if all(x in y for x, y in zip(row_series.values,
                                              self.raw_implicant))
            else 0.0, axis=1).sum() / tmp_data.apply(
            lambda row_series:
            1.0 if all(x in y for x, y in zip(row_series.values,
                                              self.raw_implicant))
            else 0.0, axis=1).sum()
        return self.incl_score
