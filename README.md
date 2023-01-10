# stPaths_assessment
This repository contains the code for the algorithms of [[1]](#1). 
More specifically, we have the code for three different algorithm for assessing the number of st-paths in an undirected graph. The algorithms will take as input a graph, a source node s, a target node t, and a threshold z. 
The algorithms will assess whether the number of paths in the given graph from node s to node t are at least z, or not, and answer accordingly.
The repository contains the code for three different algorithms:
- `Johnson-assess.cpp`, which is Johnson's algorithm for listing paths [[2]](#2), adapted to stop whenever the paths found are at least z
- `classicBCC-assess.cpp`, a simpler implementation the optimal enumeration algorithm of [[3]](#3), again adapted to stop when the paths found are at least z
- `paths-assessment.cpp`, the true assessment algorithm proposed in [[1]](#1).

Each `.cpp` file is compilable on its own. 
After regular compilation such as `g++ -O3 whichever-assess.cpp -o whichever`, the output executable will take four mandatory parameters as input; in order:
1. the input graph, provided as the directory of file representing a graph in nde format.
2. two integers, corresponding to the source node and the target node, in this order. The two nodes must belong to the same connected component, otherwise an error is raised.
3. the threshold z for which we wish to assess

The first two algorithms take two more possible optional parameters; in order:
- MAX_TIME, which sets a timeout in milliseconds for the given algorithm;
- output_filename, a file for reporting the output parameters.

Algorithm `paths-assessment.cpp` has one more optional parameter, representing the _strategy_ for source choice (see [[1]](#1) for more details). 
More specifically, the optional parameters for this algorithm are, in order, 
- MAX_TIME, which sets a timeout in milliseconds for the given algorithm;
- a character 'c' representing the strategy that the algorithm will adopt;
- output_filename, a file for reporting the output parameters.

The possible strategies for `paths-assessment.cpp` are five:
1. random, represented by character 'r'
2. smallest, represented by character 'm'
3. biggest, represented by character 'x'
4. maximum ancestor bound, represented by character 'b'
5. LIFO, represented by character 'l'

For instance, if we want to use algorithm `paths-assessment.cpp` to assess the number of paths from node 5 to node 782 of graph `inputgraph.nde` with timeout of 10 minutes, strategy given by maximum ancestor bound, and output file `outputlog.txt` is given by:  
```
g++ -O3 paths-assessment.cpp -o pathassess
./pathassess ./inputgraph.nde 5 782 600 'b' outputlog.txt 
```

When terminating, the algorithms will report the following parameters:  elapsed time, number of paths assessed, number of calles performed, maximum memory employed (MB)

It is useful to employ algorithm `preprocess.cpp` to cleanup the graph before assessment. This algorithm takes as input a graph file `input_graph.nde`, and performs the following:
1. tests pairs of nodes at random until it finds two of them that reach each other: s and t
2. removes from the graph the nodes that are not reachable from either s or t
3. outputs the cleaned graph in a new file, called `input_graph-s-t.nde`.


The files `create-circle.cpp`, `create-diamond.cpp`, and `create-ladder.cpp` can be used to create graphs of specific shapes (respectively, a cycle with several chords, a diamond graph as described in [[2]](#2), and a  ladder-shaped graph). They take in input the number of nodes and edges of the graph, and create an `.nde` file containing the graph. 

Algorithm `find-st.cpp` takes as input a graph in `.nde` format and it randomly tests pairs of nodes to find valid sources and targets (i.e. there is at least one path connecting the source to the target). It can be normally compiled as `g++ -O3 find-st.cpp -o executable_name`; the resulting executable takes just one input, the FILEDIRECTORY of the desired input graph. Thus, to execute it:
```
./executable_name ./dataset/inputgraph.nde
```
When terminating, it will output the indices of two nodes that can be used as source and target for the other algorithms.


## References
<a id="1">[1]</a> 
G. Punzi, A. Conte, R. Grossi, and A. Marino.
_An Efficient Algorithm for Assessing the Number of st-Paths in Large Graphs._
SIAM International Conference on Data Mining (SDM23), 2023, to appear.

<a id = "2">[2]</a> 
D.B. Johnson.
_Finding all the Elementary Circuits of a Directed Graph._
SIAM Journal on computing, 4(1):77–84, 1975.

<a id = "3">[3]</a> 
E. Birmelé, R. Ferreira, R. Grossi, A. Marino,
N. Pisanti, R. Rizzi, and G. Sacomoto. 
_Optimal listing of cycles and st-paths in undirected graphs._ 
ACM-SIAM SODA, pages 1884–1896. SIAM, 2013.
