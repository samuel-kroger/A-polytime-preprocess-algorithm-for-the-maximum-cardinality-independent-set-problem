# A polytime preprocess algorithm for the maximum cardinality independent set problem

This repository contains the code, data, and results for "A polytime preprocess algorithm for the maximum cardinality independent set problem" by Samuel Kroger, Hamidreza Validi, and Illya V. Hicks.

An independent set $I \subseteq V$ of a graph $G=(V,E)$ is a set such that the subgraph induced by $I$ (i.e., $G[I]$) has no edges. The maximum cardinality independent set problem aims to find a largest independent set in a graph. The maximum cardinality independent set problem is NP-hard. We propose a polytime algorithm that can fix binary decision variables of a classical mixed integer programming formulation of the problem. The preprocessing algorithm is illustrated for the karate graph below. The proposed algorithm runs in two iterations; the red nodes denote a maximum independent set of simplicial nodes, which we can fix to one; the yellow nodes are the neighbors of the red nodes that we can be fixed to zero, and the blue nodes are the remaining nodes.

![Figure 1](images/karate.jpg "The 4-core of the karate graph")



EXPLANATION OF HOW TO RUN CODE AND WHAT IS IN THE DIRECTORIES GOES HERE


## How to run the code?

```
C:\src\comp_expermient.py
```
