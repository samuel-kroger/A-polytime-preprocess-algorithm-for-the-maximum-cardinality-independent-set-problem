# A polynomial time preprocess algorithm for solving the maximum stable set problem

In this repository you will find the code, data, and results for "[placeholder](placeholder)" by Samuel Kroger, Hamidreza Validi, and Illya V. Hicks.

A stable set $I$ of a a graph is a set such that the subgraph induced by I, $G(I)$ has no edges. The maximum stable set problem is finding a largest cardinality stable set in a graph. The maximum stable set problem is NP hard. In this work we propose a polynomial time algorithm which can fix variables for the maximum stable set problem. Below you can see our preprocessing algorithm. The proposed algorithm runs in two iterations; the red nodes denote a maximum independent set of simplicial nodes which we can fix, the yellow nodes are the neighbors of the set of the red nodes which we can also fix, and the blue nodes are the renaming nodes.

![Figure 1](images/karate.jpg "The 4-core of the karate graph")



EXPLANATION OF HOW TO RUN CODE AND WHAT IS IN THE DIRECTORIES GOES HERE


## How to run the code?

```
C:\src\comp_expermient.py
```