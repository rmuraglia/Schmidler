``Schimdler'' repository contains codebase for various projects undertaken while working in the research group of Scott Schmidler at Duke University.

Work focuses on path optimization in free energy calculations. Determining minimum variance paths for free energy estimation with bridge sampling methods, and devising new estimation methods compatible with nonequilibrium sampling and edge-wise path searching are the primary goals.

Code is primarily written in R, with some bash scripts to facilitate batch running. Some code for earlier projects are trimmed down versions from their original copies that were developed outside of a version control setting. There may potentially be some bugs, but on the whole, this code represents, at minimum, a precise blueprint for our studies.

The contents of each directory are detailed below: 

* _Thesis_: TeX files for Master's thesis. Includes figures and defense slides.
* _Exhaustive-search_: Code for chapter 3 (first results chapter) of thesis. Mix of Dijkstra search and other dynamic programming techniques for optimal path finding.
* _SMC_: Code for chapter 4 and first figure of chapter 5. Details use of Sequential Monte Carlo methods for ratio estimation. Majority of work done with AIS, implementation details for seqBAR and pCrooks.
* _Q-learning_: Code for chapter 5. Q-learning for determining minimum variance path for ratio estimation with sequential sampling (pCrooks). 

