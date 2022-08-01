# MOTIF-DISCOVERY_GIBBS-SAMPLER-using-python
Gibbs sampling (also called alternating conditional sampling) is a Markov Chain Monte Carlo algorithm for high-dimensional data.
GibbsSampler is a motif finding algorithm that finds one common motif and returns a list of bestMotifs containing the closest motif match from each string in dna.
# <strong>Problem Statement</strong>
Given: Integers k, t, and N, followed by a collection of strings Dna.<br>
Return: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts.<br>
We have to code Gibbs Sampler for the purpose of motif discovery. Here,<br>
    Dna -- A collection of DNA strings that are of the same length.<br>
    "t" -- Is an integer indicating how many times to read the genetic algorithm.<br>
    "k" -- An integer indicating the motif length being searched for.<br>
    "N" -- The number of iterations before returning the best motif.<br>
    
# Algorithm
GibbsSampler (Dna, k , t , N)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;randomly select k−Mers Motifs = ( Motif 1 , . . . , Motift ) in each string from Dna<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;BestMotifs ←Motifs<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for j ←1 to N<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;i ← Random ( t )<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Profile ←profile matrix constructed from all strings in Motifs except for Motifi<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Motifi ←Profile − randomly generated k-Mer in the i-th sequence<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;if Score ( Motifs ) < Score ( BestMotifs )<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;BestMotifs ←Motifs<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;return BestMotifs <br>

