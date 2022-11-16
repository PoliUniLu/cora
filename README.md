# CORA

CORA is a Python library for **Combinational Regularity Analysis** (CORA). 

## Description

CORA is a member of the family of Configurational Comparative Methods (CCMs; Thiem *et al*. 2016). It is thus closely related to Qualitative Comparative Analysis (QCA; Ragin 1987) and Coincidence Analysis (CNA; Baumgartner 2009; Baumgartner and Ambühl 2020). CCMs seek to detect so-called "INUS" structures in data (Thiem 2017, 2022). Such structures represent patterns of causal relations modelled in the language of propositional logic (Mackie 1965; Psillos 2009). 

In contrast to QCA and CNA, however, CORA is centrally inspired by the field of switching circuit analysis, a subfield of electrical engineering. The reason is that *propositional logic* - the language of INUS causation - and *switching algebra* - the language of switching circuit analysis - are two conceptually distinct yet operationally equivalent branches of the same underlying Boolean algebra (Lewin and Protheroe 1992). As such, CORA is currently the only CCM that can analyze multi-output structures, that is, INUS structures that feature complex conjunctive effects. In addition, CORA offers a data-mining approach to solution building that reduces model ambiguities by keeping the number of required variables for finding a solution at a minimum. Lastly, CORA includes a visualization module called LOGIGRAM, with which logic diagrams can be produced from any (system of) Boolean function(s) in disjunctive normal form.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install CORA.

```bash
pip install CORA
```
## Google Colab

To open CORA in Google Colab, click the button below:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/gist/ZuzanaSebb/490501c7573490ebfbb841e850b573f9/-cora_1-0-3.ipynb)

## Citation Info

When using CORA, please cite as follows:

Sebechlebská, Zuzana, Lusine Mkrtchyan and Alrik Thiem. 2022. CORA: A Python package for Combinational Regularity Analysis, Version *current version number*. Available from: https://github.com/PoliUniLu/cora.


## Copyright

CORA is licensed under a CC BY-NC-SA 4.0 license. This license allows reusers to distribute, remix, adapt, and build upon the material in any medium or format for noncommercial purposes only, and only so long as attribution is given to the creator. If you remix, adapt, or build upon the material, you must license the modified material under identical terms. 

[![Open In License](https://mirrors.creativecommons.org/presskit/buttons/88x31/svg/by-nc-sa.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)


## References

* Baumgartner, Michael. 2009. "Inferring Causal Complexity." *Sociological Methods & Research* **38** (1):71-101. https://doi.org/10.1177/0049124109339.
* Baumgartner, Michael, and Mathias Ambühl. 2020. "Causal modeling with multi-value and fuzzy-set Coincidence Analysis." *Political Science Research and Methods* **8** (3):526-42. https://doi.org/10.1017/psrm.2018.45.
* Lewin, Douglas, and David Protheroe. 1992. *Design of Logic Systems*. 2nd ed. London: Chapman & Hall.
* Mackie, John L. 1965. "Causes and Conditions." *American Philosophical Quarterly* **2** (4):245-64. https://www.jstor.org/stable/20009173.
* Psillos, Stathis. 2009. "Regularity Theories." In *The Oxford Handbook of Causation*, ed. H. Beebee, C. Hitchcock and P. Menzies. Oxford: Oxford University Press, pp.131-157.
* Ragin, Charles C. 1987. *The Comparative Method: Moving beyond Qualitative and Quantitative Strategies*. Berkeley: University of California Press.
* Thiem, Alrik. 2017. "Conducting Configurational Comparative Research with Qualitative Comparative Analysis: A Hands-On Tutorial for Applied Evaluation Scholars and Practitioners." *American Journal of Evaluation* **38** (3):420-33. https://doi.org/10.1177/109821401667.
* Thiem, Alrik. 2022. "Qualitative Comparative Analysis (QCA)." In *Handbook of Research Methods in International Relations*, ed. R. J. Huddleston, T. Jamieson and P. James. Cheltenham: Edward Elgar, pp.607-28. https://doi.org/10.4337/9781839101014.00044.
* Thiem, Alrik, Michael Baumgartner, and Damien Bol. 2016. "Still Lost in Translation! A Correction of Three Misunderstandings Between Configurational Comparativists and Regressional Analysts." *Comparative Political Studies* **49** (6):742-74. https://doi.org/10.1177/00104140145658.
