import itertools
import pandas as pd
import cora


class TupleResult:

    def __init__(
        self, combination, out_len, irrendudant_systems, inc_score, cov_score, score
    ):

        self.combination = combination
        if (
            out_len == 1
            and len(irrendudant_systems) == 1
            and (str(irrendudant_systems[0].system[0].implicant) == "1")
        ):
            self.nr_irr_systems = 0
            self.inc_score = 0
            self.cov_score = 0
            self.score = score = 0
        else:
            self.nr_irr_systems = len(irrendudant_systems)
            self.inc_score = inc_score
            self.cov_score = cov_score
            self.score = score


def data_mining(
    data,
    output_labels,
    len_of_tuple,
    input_labels=None,
    case_col=None,
    n_cut=1,
    inc_score1=1,
    inc_score2=None,
    U=None,
    algorithm="ON-DC",
    automatic=False,
):
    """
    Function performing data mining on a subset of input variables.

    Parameters
    ----------
    data : dataframe

    output_labels : an array of strings
    The names of the outcome columns from the data frame.

    len_of_tuple : int
    Number indicating how many input variables from the original data
    are used in the computation.

    case_col : string
    The name of the column from the data frame containing the case ids.

    input_labels : an array of strings
    The names of the input columns from the data frame.

    n_cut : int
    The minimum number of cases under which a truth table row is declared as a
    don't care.

    inc_score1 : float
    The minimum sufficiency inclusion score for an output function value of "1".

    inc_score2 : float
    The maximum sufficiency inclusion score for an output function value of "0".

    U : int
    The U number is either 0 or 1.

    algorithm : string
    The name of the optimization algorithm

    automatic : boolean
    If true the function will return the first non-zero solution approached by
    the data-mining search.

    Returns
    -------
    result : dataframe

    Example
    -------
    >>> data = pd.DataFrame([[1,0,1,0],
    ...                     [1,1,1,1],
    ...                     [1,0,0,0],
    ...                     [0,1,0,1]],
    ...                     columns=["A","B","C","O"])
    >>> data_mining(data,["O"],2))

        Combination  Nr_of_systems  Inc_score  Cov_score  Score
    0      [A, B]              1        1.0        1.0    1.0
    1      [A, C]              1        1.0        0.5    0.5
    2      [B, C]              1        1.0        1.0    1.0
    """
    res = []

    if input_labels is None:
        input_labels = [
            x for x in data.columns if (x not in output_labels and x != case_col)
        ]
    for i, comb in enumerate(
        itertools.combinations([x for x in input_labels], len_of_tuple)
    ):
        cols = list(comb)
        data_object = cora.OptimizationContext(
            data,
            output_labels,
            cols,
            case_col=case_col,
            n_cut=n_cut,
            inc_score1=inc_score1,
            inc_score2=inc_score2,
            U=U,
            algorithm=algorithm,
        )

        if len(output_labels) == 1:
            ir_sys = data_object.get_irredundant_sums()
        else:
            ir_sys = data_object.get_irredundant_systems()

        if len(ir_sys) > 0:
            inc_score = round(max(x.inclusion_score() for x in ir_sys), 3)
            cov_score = round(max(x.coverage_score() for x in ir_sys), 3)
            score = round(
                max(x.inclusion_score() * x.coverage_score() for x in ir_sys), 3
            )
            tr = TupleResult(
                cols, len(output_labels), ir_sys, inc_score, cov_score, score
            )
            res.append(tr)

        else:
            res.append(TupleResult(cols, len(output_labels), ir_sys, 0, 0, 0))

    rows = [
        (x.combination, x.nr_irr_systems, x.inc_score, x.cov_score, x.score)
        for x in res
    ]
    result = pd.DataFrame(
        rows,
        columns=["Combination", "Nr_of_systems", "Inc_score", "Cov_score", "Score"],
    )

    if (
        not automatic
        or sum(result["Nr_of_systems"]) != 0
        or len_of_tuple >= len(input_labels)
    ):
        return result

    else:
        return data_mining(
            data,
            output_labels,
            len_of_tuple + 1,
            input_labels,
            case_col,
            n_cut,
            inc_score1,
            inc_score2,
            U,
            algorithm,
            automatic,
        )


if __name__ == "__main__":
    df = pd.DataFrame(
        [[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 0], [1, 0, 1, 1]],
        columns=["A", "B", "C", "Z"],
    )
    print(data_mining(df, ["Z"], 1, automatic=True))
