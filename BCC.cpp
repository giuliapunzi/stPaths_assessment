#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>

using namespace std;

int N; // size of the graph/biggest component; used for global vectors
int maxID = 0; // maximum ID of a BCC

typedef struct BCC {
    vector<int> sources; // vector of paths' sources contained in the BCC
    int target; // target of the BCC = only articulation point leading to t 
    long personalBound; // personal bound of current BCC = |edges| - |nodes| + 1
    long prodAncestors; // product of the personal bounds of all BCCs in the bead string from current to root
    int myID; // personal ID of the BCC
    bool isLeaf; // true iff the BCC is a leaf of the block-cut tree

    vector<int> nodes; // vector of node indices forming the current BCC
    unordered_map<int, vector<int>> edges; // unordered map assigning to every node its incident edges
};


// these are global vectors 
vector<bool> visited;
vector<int> parent;
vector<int> tree_stack;
vector<int> disc;
vector<int> low;
int visit_time;
vector<bool> is_art_point;

void findBCCs(int u, BCC &B, vector<BCC> &BCC_vector, vector<int> &source_neighbors)
{
    tree_stack.push_back(u);
    // Count of children in DFS Tree (for root)
    int children = 0;
 
    // Mark the current node as visited
    visited[u] = true;
 
    // Initialize discovery time and lowpoint value
    disc[u] = low[u] = ++visit_time;

    // Go through all neighbors of u
    for (auto v : B.edges.at(u)) {
        // if v is not visited yet, then make it a child of u in DFS tree and recur for it
        if (!visited[v]) {
            parent[v] = u;
            children++;
            findBCCs(v, B, BCC_vector, source_neighbors);

            // if we are the root, we just closed a BCC
            // if(parent[u] == -1){
            //     BCC current_BCC;
            //     current_BCC.target = u; // u is the target of the current BCC

            //     // for the current BCC we need to unstack until v
            //     int x = tree_stack.back();
            //     while(x != v){ // If we are in the "right" BCC add it to the vector   
            //         current_BCC.nodes.push_back(x);
            //         tree_stack.pop_back();
            //         x = tree_stack.back();
            //     }
            //     tree_stack.pop_back(); 
            //     current_BCC.nodes.push_back(v);

            //     current_BCC.isLeaf = false;
            //     current_BCC.sources = {};
            //     current_BCC.myID = ++ID;
            //     current_BCC.prodAncestors = -1;

            //     // we fill edges right away; other structure is filled at the end
            //     // for each node x, iterate through incident edges for it in B and add edges if the other endpoint is also in currentBCC
            //     for(auto x: current_BCC.nodes) {
            //         current_BCC.edges.insert({x, {}}); // start with empty vector for all nodes
            //         for (auto y: B.edges.at(x))
            //         {
            //             if(find(current_BCC.nodes.begin(), current_BCC.nodes.end(), y) != current_BCC.nodes.end()){
            //                 current_BCC.edges[x].push_back(y); // add y to edges' vector of x
            //             }
            //         }
            //     }

            //     // add the BCC to the final vector
            //     BCC_vector.push_back(current_BCC);
            // }

            // Check if the subtree rooted with v has
            // a connection to one of the ancestors of u
            low[u] = min(low[u], low[v]);

            // if u is not root and low value of one of its child is more than discovery value of u, 
            // OR if u is the root (returning from a child of the root identifies a BCC), then we close a BCC
            if (parent[u] != -1 && low[v] >= disc[u]){ // here is where I close my articulation point
                BCC current_BCC;
                current_BCC.target = u; // u is the target of the current BCC
                is_art_point[u] = true;

                // for the current BCC we need to unstack until v
                int x = tree_stack.back();
                while(x != v){ // If we are in the "right" BCC add it to the vector   
                    current_BCC.nodes.push_back(x);
                    // check if x is an art point; if it is, add it to sources
                    if(is_art_point[x])
                        current_BCC.sources.push_back(x);

                    // also check if x is a neighbor of the source
                    if(find(source_neighbors.begin(), source_neighbors.end(), x) != source_neighbors.end())
                        current_BCC.sources.push_back(x);

                    tree_stack.pop_back();
                    x = tree_stack.back();
                }
                tree_stack.pop_back(); 
                if(is_art_point[v])
                    current_BCC.sources.push_back(v);

                // also check if x is a neighbor of the source
                if(find(source_neighbors.begin(), source_neighbors.end(), v) != source_neighbors.end())
                    current_BCC.sources.push_back(v);

                current_BCC.nodes.push_back(v); // always finish at v

                current_BCC.isLeaf = false;
                current_BCC.sources = {};
                current_BCC.myID = ++maxID;
                current_BCC.prodAncestors = -1;

                // for each node x, iterate through incident edges for it in B and add edges if the other endpoint is also in currentBCC
                for(auto x: current_BCC.nodes) {
                    current_BCC.edges.insert({x, {}}); // start with empty vector for all nodes
                    for (auto y: B.edges.at(x))
                    {
                        if(find(current_BCC.nodes.begin(), current_BCC.nodes.end(), y) != current_BCC.nodes.end()){
                            current_BCC.edges[x].push_back(y); // add y to edges' vector of x
                        }
                    }
                }

                // add the BCC to the final vector
                BCC_vector.push_back(current_BCC);

            }
        }

        // Update low value of u for parent function calls.
        else if (v != parent[u])
            low[u] = min(low[u], disc[v]);
    }
    
    return;
}



// INPUT: a BCC B and one of its sources s
// OUTPUT: a vector of non trivial BCCs which form the block-cut tree of B after the removal of s
vector<BCC> expand(BCC B, int s){
    // check if s is a source of to B
    if(find(B.sources.begin(), B.sources.end(), s) == B.sources.end())
        throw invalid_argument("Source provided is not a source for the BCC");

    // remove s from nodes, incident edges of s from edges
    B.nodes.erase(find(B.nodes.begin(), B.nodes.end(), s));
    vector<int> source_neighbors = B.edges[s]; // save neighbors of s for later
    B.edges.erase(B.edges.find(s));

    // initialize vectors (maybe make them global)
    tree_stack.erase(tree_stack.begin(), tree_stack.end());
    visit_time = 0;

    for(int i = 0; i < N; i++){
        visited[i] = false;
        parent[i] = -2;
        is_art_point[i] = false;
    }
    
    parent[B.target] = -1;

    // start visit from target of BCC (we first close BCCs that contain neighbors of s)
    vector<BCC> decomposed_BCC = {};
    findBCCs(B.target, B, decomposed_BCC, source_neighbors);

    // after call for art points, finish unstacking and form last new BCC
    BCC current_BCC;
    current_BCC.target = B.target; // B.target is the target of the current BCC
    is_art_point[B.target] = true;

    // for the current BCC we need to unstack until the stack is empty
    int x = tree_stack.back();
    while(!tree_stack.empty()){ 
        current_BCC.nodes.push_back(x);
        if(is_art_point[x])
            current_BCC.sources.push_back(x);

        // also check if x is a neighbor of the source
        if(find(source_neighbors.begin(), source_neighbors.end(), x) != source_neighbors.end())
            current_BCC.sources.push_back(x);
        
        tree_stack.pop_back();
        x = tree_stack.back();
    }
    tree_stack.pop_back(); 
    if(is_art_point[x])
        current_BCC.sources.push_back(x);

    // also check if x is a neighbor of the source
    if(find(source_neighbors.begin(), source_neighbors.end(), x) != source_neighbors.end())
        current_BCC.sources.push_back(x);

    current_BCC.nodes.push_back(x);

    current_BCC.isLeaf = false;
    current_BCC.sources = {};
    current_BCC.myID = ++maxID;
    current_BCC.prodAncestors = -1;

    // fill edges right away: for each node x, iterate through incident edges for it in B and add edges if the other endpoint is also in currentBCC
    for(auto x: current_BCC.nodes) {
        current_BCC.edges.insert({x, {}}); // start with empty vector for all nodes
        for (auto y: B.edges.at(x))
        {
            if(find(current_BCC.nodes.begin(), current_BCC.nodes.end(), y) != current_BCC.nodes.end()){
                current_BCC.edges[x].push_back(y); // add y to edges' vector of x
            }
        }
    }

    // add the BCC to the final vector
    decomposed_BCC.push_back(current_BCC);

    // need to set up personal bounds, isLeaf parameter by iterating through all new components
    for (auto comp : decomposed_BCC)
    {   
        comp.personalBound = comp.edges.size() - comp.nodes.size() +1; // personal bound given by cyclomatic lemma
        
        // it is a leaf if the sources vector has size one, and the only node belongs to the neighbors of s
        if(comp.sources.size() == 1 && find(source_neighbors.begin(), source_neighbors.end(), comp.sources[0])!= source_neighbors.end())
            comp.isLeaf = true;
    }
    

    return decomposed_BCC;
}
