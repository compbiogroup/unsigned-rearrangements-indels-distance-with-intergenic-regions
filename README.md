### Supplemental Material: Incorporating Intergenic Regions into Reversal and Transposition Distances with Indels

The file *Supplemental Material.pdf* presents examples of the main concepts used in the algorithms proposed in the paper and experimental results using simulated data.

The folder *code/Instances* has a script to generate a set of random instances. Furthermore, it has the datasets used in the experiments reported in the supplemental material of the paper.

The approximation algorithms described in the paper are in the folder *code/algorithms*.

Important: The results presented in the *Supplemental Material.pdf* for Reversals, Transpositions, and Indels use the 6-approximation algorithm, which is the same algorithm as the 4-approximation for Reversals and Indels. We developed a new algorithm, which we added to the *code* folder, but it is not the algorithm presented in the paper. To reproduce the results from the paper, change the file *main.py* to use the Reversals and Indels algorithm for the model with Reversals, Transpositions, and Indels. 
