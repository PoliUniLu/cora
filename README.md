# CORA
[![Tests](https://github.com/PoliUniLu/cora/workflows/Cora/badge.svg)](https://github.com/PoliUniLu/cora/actions?workflow=Cora)
[![codecov](https://codecov.io/github/PoliUniLu/cora/branch/master/graph/badge.svg?token=36V7QBSJI3)](https://codecov.io/github/PoliUniLu/cora)
[![DOI](https://zenodo.org/badge/289114682.svg)](https://zenodo.org/badge/latestdoi/289114682)



CORA is a Python library for **Combinational Regularity Analysis** (CORA). 

## Description

Combinational Regularity Analysis (CORA) is a member of the family of Configurational Comparative Methods (CCMs; Thiem *et al*. 2022). It is thus closely related to Qualitative Comparative Analysis (QCA; Ragin 1987) and Coincidence Analysis (CNA; Baumgartner 2009; Baumgartner and Ambühl 2020). Modern CCMs seek to detect INUS structures in data (Thiem 2017, 2022). Such structures are elaborate cause-effect relations that can be represented in the Boolean language of propositional logic (Baumgartner 2008; Mackie 1965; Psillos 2009). Most importantly, these relations are marked by causal conjunctivity (e.g., *a* **and** **not** *b* **and** *c* **and** $\cdots$) and causal disjunctivity (e.g., *d* **or** *e* **or** *f* **or** $\cdots$). For this reason, CCMs differ fundamentally from most other empirical research methods (Thiem *et al*. 2016).

In contrast to QCA and CNA, however, CORA has been inspired by switching circuit analysis, a subfield of electrical engineering. INUS structures and switching circuits have much in common because propositional logic - the language of INUS causation - and switching algebra - the language of switching circuit analysis - are operationally equivalent branches of the same underlying Boolean algebra (Lewin and Protheroe 1992). It is therefore no coincidence that one of the first systematic algorithms for Boolean optimization - the Quine-McCluskey algorithm (McCluskey 1956; Quine 1955) - had been co-developed by an analytical philosopher (Willard Van Orman Quine) and an electrical engineer (Edward J. McCluskey).

Most importantly, CORA is currently the only CCM able to analyze INUS structures that simultaneously feature simple as well as complex effects (e.g., *y* **and** **not** *z*, **not** *y* **and** *z*, *y* **and** *z*). CORA can process such structures even in multi-value form (Mkrtchyan *et al*. 2023). In addition, CORA offers a configurational version of Occam's Razor: a data-mining approach to solution building that reduces model ambiguities by keeping the number of required variables for finding a solution at a minimum. Lastly, CORA includes a lean yet powerful visualization module called LOGIGRAM, with which two-level logic diagrams can be produced from any (system of) Boolean or multi-value function(s) in disjunctive normal form. Logic diagrams considerably outperfrom Venn diagrams, which are often used in QCA, when it comes to the representation and interpretability of INUS structures (Thiem *et al*. 2023).

## Installation

CORA requires Python>3.7 and uses Poetry for dependendecy management. Use the package manager [pip](https://pip.pypa.io/en/stable/) to install CORA, including all dependencies.

```bash
pip install git+https://github.com/PoliUniLu/cora.git
```

It is recommended to install the package into a dedicated virtual environment.

## Google Colab

To open CORA with a graphical interface in Google Colab, click the button below:

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/gist/ZuzanaSebb/5a47670cf0e19126ded51127e78270e8/-cora_1-0-3.ipynb)

## Usage
The main features of the package are part of the `OptimizationContext` class, including functions:
- `get_prime_implicants`,
- `prime_implicant_chart`,
- `get_irredundant_systems`,
- `get_irredundant_solutions`.

***Note:***
Use the `help` function to access the documentation of CORA.

**Example:**
```python
import pandas as pd
from cora import OptimizationContext

df = pd.DataFrame([[1,1,0,1],
                   [0,0,1,1],
                   [1,0,1,0],
                   [0,1,0,1]], columns=["A","B","C","OUT"])

context = OptimizationContext(data = df, output_labels = ["OUT"])
PIs = context.get_prime_implicants() # result: {B, c, #a}; essential prime implicants marked by hashtags
irredundant_solutions = context.get_irredundant_sums() # result: [M1: #a + B, M2: #a + c]
```
Configurational data-mining is another feature. It analyzes all n-tuples of input combinations to search for feasible tuples of solution-generating inputs. In essence, this feature thus provides a configurational version of Occam's Razor (Feldman 2016). 

**Example:**
```python
import pandas as pd
import cora

df = pd.DataFrame([[1,2,0,1,1],
                   [1,1,1,0,1],
                   [0,2,1,0,0],
                   [0,2,2,0,1],], columns=["A","B","C","D","OUT"])
result = cora.data_mining(df, ["OUT"], len_of_tuple = 2, inc_score1 = 0.5, n_cut = 1)
result # print(result.to_markdown())

|    | Combination   |   Nr_of_systems |   Inc_score |   Cov_score |   Score |
|---:|:--------------|----------------:|------------:|------------:|--------:|
|  0 | ['A', 'B']    |               1 |        0.75 |           1 |    0.75 |
|  1 | ['A', 'C']    |               1 |        1    |           1 |    1    |
|  2 | ['A', 'D']    |               1 |        0.75 |           1 |    0.75 |
|  3 | ['B', 'C']    |               1 |        1    |           1 |    1    |
|  4 | ['B', 'D']    |               1 |        0.75 |           1 |    0.75 |
|  5 | ['C', 'D']    |               1 |        0.75 |           1 |    0.75 |

```
To access more examples, see the `/examples` folder or follow 
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/PoliUniLu/cora/blob/master/examples/cora_examples_notebook.ipynb)



## Citation Info

When using CORA (method and software), please cite it as follows:

**Method:**

Thiem, Alrik, Lusine Mkrtchyan, and Zuzana Sebechlebská. 2022. "Combinational Regularity Analysis (CORA) - A New Method for Uncovering Complex Causation in Medical and Health Research." *BMC Medical Research Methodology*  **22** (1):333. [Link](http://dx.doi.org/10.1186/s12874-022-01800-9)

**Software:**

Sebechlebská, Zuzana, Lusine Mkrtchyan and Alrik Thiem. 2022. CORA: A Python package for Combinational Regularity Analysis, Version <*current version number*>. Available from: https://github.com/PoliUniLu/cora.


## Copyright

CORA is licensed under GNU GPLv3. 

## Contributions 

We highly welcome contributions from the community. Feedback, bug reports, and feature requests should be placed as an [Issue](https://github.com/PoliUniLu/cora/issues) on GitHub.

### Pull requests
To set up a development environment, use [Poetry](https://python-poetry.org/).
```console
pip install poetry
poetry install
```
Test the code by running
```console
poetry run pytest
```
Pull requests are welcome. Note that, although the current codebase does not have an entirely consistent code style, the new code should be PEP-8 compliant.

## References

* Baumgartner, Michael. 2008. "Regularity Theories Reassessed." *Philosophia* **36** (3):327-54. [Link](https://doi.org/10.1177/0049124109339)
* Baumgartner, Michael. 2009. "Inferring Causal Complexity." *Sociological Methods & Research* **38** (1):71-101. [Link](https://doi.org/10.1177/0049124109339)
* Baumgartner, Michael, and Mathias Ambühl. 2020. "Causal modeling with multi-value and fuzzy-set Coincidence Analysis." *Political Science Research and Methods* **8** (3):526-42. [Link](https://doi.org/10.1017/psrm.2018.45)
* Feldman, Jacob. 2016. "The Simplicity Principle in Perception and Cognition." *Wiley Interdisciplinary Reviews: Cognitive Science* **7** (5):330-40. [Link](https://doi.org/10.1002/wcs.1406)
* Lewin, Douglas, and David Protheroe. 1992. *Design of Logic Systems*. 2nd ed. London: Chapman & Hall.
* Mackie, John L. 1965. "Causes and Conditions." *American Philosophical Quarterly* **2** (4):245-64. [Link](https://www.jstor.org/stable/20009173)
* McCluskey, Edward J. 1956. "Minimization of Boolean Functions." *Bell Systems Technical Journal* **35** (6):1417-44. [Link](http://onlinelibrary.wiley.com/doi/10.1002/j.1538-7305.1956.tb03835.x/abstract)
* Mkrtchyan, Lusine, Alrik Thiem, and Zuzana Sebechlebská. 2023. "Re-Establishing a Lost Connection: Multi-Value Logic in Causal Data Analysis in Social Science Disciplines." *IEEE Access* **11**:10471-82. [Link](https://doi.org/10.1109/ACCESS.2023.3240094)
* Psillos, Stathis. 2009. "Regularity Theories." In *The Oxford Handbook of Causation*, ed. H. Beebee, C. Hitchcock and P. Menzies. Oxford: Oxford University Press, pp.131-157. [Link](https://doi.org/10.1093/oxfordhb/9780199279739.003.0008)
* Quine, Willard V. O. 1955. "A Way to Simplify Truth Functions." *American Mathematical Monthly* **62** (9):627-31. [Link](http://www.jstor.org/stable/2307285)
* Ragin, Charles C. 1987. *The Comparative Method: Moving beyond Qualitative and Quantitative Strategies*. Berkeley: University of California Press. [Link](https://www.jstor.org/stable/10.1525/j.ctt1pnx57)
* Thiem, Alrik. 2017. "Conducting Configurational Comparative Research with Qualitative Comparative Analysis: A Hands-On Tutorial for Applied Evaluation Scholars and Practitioners." *American Journal of Evaluation* **38** (3):420-33. [Link](https://doi.org/10.1177/109821401667)
* Thiem, Alrik. 2022. "Qualitative Comparative Analysis (QCA)." In *Handbook of Research Methods in International Relations*, ed. R. J. Huddleston, T. Jamieson and P. James. Cheltenham: Edward Elgar, pp.607-28. [Link](https://doi.org/10.4337/9781839101014.00044)
* Thiem, Alrik, Lusine Mkrtchyan, and Zuzana Sebechlebská. 2022. "Combinational Regularity Analysis (CORA) - A New Method for Uncovering Complex Causation in Medical and Health Research." *BMC Medical Research Methodology*  **22** (1):333. [Link](http://dx.doi.org/10.1186/s12874-022-01800-9)
* Thiem, Alrik, Zuzana Sebechlebská, and Lusine Mkrtchyan. 2023. "Logic Diagrams: A Visual Tool with Untapped Potential." *Nature Reviews Methods Primers* **3** (1):21. [Link](https://rdcu.be/c7dys)
* Thiem, Alrik, Michael Baumgartner, and Damien Bol. 2016. "Still Lost in Translation! A Correction of Three Misunderstandings Between Configurational Comparativists and Regressional Analysts." *Comparative Political Studies* **49** (6):742-74. [Link](https://doi.org/10.1177/00104140145658)
