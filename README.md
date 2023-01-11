# CORA

CORA is a Python library for **Combinational Regularity Analysis** (CORA). 

## Description

CORA is a member of the family of Configurational Comparative Methods (CCMs; Thiem *et al*. 2022). It is thus closely related to Qualitative Comparative Analysis (QCA; Ragin 1987) and Coincidence Analysis (CNA; Baumgartner 2009; Baumgartner and Ambühl 2020). Modern CCMs seek to detect INUS structures in data (Thiem 2017, 2022). Such structures are specific cause-effect relations that can be represented in the Boolean language of propositional logic (Baumgartner 2008; Mackie 1965; Psillos 2009). Most importantly, these relations are marked by causal conjunctivity (*a* **and** **not** *b* **and** *c* **and** $\cdots$) and causal disjunctivity (*d* **or** *e* **or** *f* **or** $\cdots$). CCMs thus differ fundamentally from most other empirical research methods (Thiem *et al*. 2016).

In contrast to QCA and CNA, however, CORA is inspired by switching circuit analysis, a subfield of electrical engineering. INUS structures and switching circuits have much in common because propositional logic - the language of INUS causation - and switching algebra - the language of switching circuit analysis - are two conceptually distinct yet operationally equivalent branches of the same underlying Boolean algebra (Lewin and Protheroe 1992). 

CORA is currently the only CCM able to analyze INUS structures that simultaneously feature simple as well as conjunctively complex effects (e.g., *y* **and** **not** *z*, **not** *y* **and** *z*, *y* **and** *z*). In addition, CORA offers a configurational version of Occam's Razor: a data-mining approach to solution building that reduces model ambiguities by keeping the number of required variables for finding a solution at a minimum. Lastly, CORA includes a lean yet powerful visualization module called LOGIGRAM, with which two-level logic diagrams can be produced from any (system of) Boolean function(s) in disjunctive normal form.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install CORA.

```bash
pip install git+https://github.com/PoliUniLu/cora.git
```
## Google Colab

To open CORA in Google Colab, click the button below:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/gist/ZuzanaSebb/d2fb79f9466af1e74b2de7047d258388/cora-1-0-5.ipynb?authuser=1)

## Citation Info

When using CORA (method and software), please cite it as follows:

**Method:**
Thiem, Alrik, Lusine Mkrtchyan, and Zuzana Sebechlebská. 2022. "Combinational Regularity Analysis (CORA) - A New Method for Uncovering Complex Causation in Medical and Health Research." *BMC Medical Research Methodology*  **22** (1):333. [Link](http://dx.doi.org/10.1186/s12874-022-01800-9)

**Software:**
Sebechlebská, Zuzana, Lusine Mkrtchyan and Alrik Thiem. 2022. CORA: A Python package for Combinational Regularity Analysis, Version <*current version number*>. Available from: https://github.com/PoliUniLu/cora.


## Copyright

CORA is licensed under a GNU GPLv3. 


## References

* Baumgartner, Michael. 2008. "Regularity Theories Reassessed." *Philosophia* **36** (3):327-54. [Link](http://dx.doi.org/10.1007/s11406-007-9114-4)
* Baumgartner, Michael. 2009. "Inferring Causal Complexity." *Sociological Methods & Research* **38** (1):71-101. [Link](https://doi.org/10.1177/0049124109339)
* Baumgartner, Michael, and Mathias Ambühl. 2020. "Causal modeling with multi-value and fuzzy-set Coincidence Analysis." *Political Science Research and Methods* **8** (3):526-42. [Link](https://doi.org/10.1017/psrm.2018.45)
* Lewin, Douglas, and David Protheroe. 1992. *Design of Logic Systems*. 2nd ed. London: Chapman & Hall.
* Mackie, John L. 1965. "Causes and Conditions." *American Philosophical Quarterly* **2** (4):245-64. [Link](https://www.jstor.org/stable/20009173)
* Psillos, Stathis. 2009. "Regularity Theories." In *The Oxford Handbook of Causation*, ed. H. Beebee, C. Hitchcock and P. Menzies. Oxford: Oxford University Press, pp.131-157. [Link](https://doi.org/10.1093/oxfordhb/9780199279739.003.0008)
* Ragin, Charles C. 1987. *The Comparative Method: Moving beyond Qualitative and Quantitative Strategies*. Berkeley: University of California Press. [Link](https://www.jstor.org/stable/10.1525/j.ctt1pnx57)
* Thiem, Alrik. 2017. "Conducting Configurational Comparative Research with Qualitative Comparative Analysis: A Hands-On Tutorial for Applied Evaluation Scholars and Practitioners." *American Journal of Evaluation* **38** (3):420-33. [Link](https://doi.org/10.1177/109821401667)
* Thiem, Alrik. 2022. "Qualitative Comparative Analysis (QCA)." In *Handbook of Research Methods in International Relations*, ed. R. J. Huddleston, T. Jamieson and P. James. Cheltenham: Edward Elgar, pp.607-28. [Link](https://doi.org/10.4337/9781839101014.00044)
* Thiem, Alrik, Lusine Mkrtchyan, and Zuzana Sebechlebská. 2022. "Combinational Regularity Analysis (CORA) - A New Method for Uncovering Complex Causation in Medical and Health Research." *BMC Medical Research Methodology*  **22** (1):333. [Link](http://dx.doi.org/10.1186/s12874-022-01800-9)
* Thiem, Alrik, Michael Baumgartner, and Damien Bol. 2016. "Still Lost in Translation! A Correction of Three Misunderstandings Between Configurational Comparativists and Regressional Analysts." *Comparative Political Studies* **49** (6):742-74. [Link](https://doi.org/10.1177/00104140145658)
