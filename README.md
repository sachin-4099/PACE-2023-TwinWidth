# HATTER: Hybrid Aprroach for Twinwidth by Taming Red (PACE 2023)

This repository contains a **heuristic solution** to the Twinwidth problem.

The program is submitted to the [PACE Challenge 2023](https://pacechallenge.org/2023/) in the **heuristic track**.

The **Directed Feedback Vertex Set** problem is defined as follows:

**Input**: An undirected graph $G = (V, E)$

**Output**: A contraction sequence that contracts $G$ to a single vertex.

**Width**: The width of the contraction sequence.

## Compilation & Usage

The whole program is contained in a single C++ file with no dependencies beyond the STL. 

You can use it as follows:

```
g++ -std=c++14 -O2 hatter.cpp -o hatter
cat input_graph.gr | ./hatter > twinwidth.txt
```
By default, the solver runs till it receives a `SIGTERM`.

## Description
