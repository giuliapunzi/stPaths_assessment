#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>
#include <chrono>
#include <cstdint>

uint64_t timeMs() {
  using namespace std::chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

// #define DEBUG true
#define DEBUG false

#define RANDOM_NODE_EXTRACTION

using namespace std;

int N; // size of the graph/biggest component; used for global vectors
int maxID = 0; // maximum ID of a BCC

struct BCC {
    vector<pair<int,int>> sources; // vector of paths' sources contained in the BCC
    int target; // target of the BCC = only articulation point leading to t 
    long personalBound; // personal bound of current BCC = |edges| - |nodes| + 1
    long prodAncestors; // product of the personal bounds of all BCCs in the bead string from current to root
    int myID; // personal ID of the BCC
    bool isLeaf; // true iff the BCC is a leaf of the block-cut tree

    vector<int> nodes; // vector of node indices forming the current BCC
    unordered_map<int, vector<int>> edges; // unordered map assigning to every node its incident edges
    int num_edges;
};


void printBCC(BCC B){
    cout << "\tGraph: ";
    for(auto x : B.nodes){
        cout <<"\t"<< x << ": [";
        for(auto y : B.edges[x]){
            cout << y << " ";
        }
        cout << "]"; 
    }
    cout << endl;
    cout << "\tPersonal bound: " << B.personalBound << endl;
    cout << "\tID: " << B.myID << endl;
    cout << "\tSources: ";
    for(auto ss : B.sources)
        cout << "(" << ss.first << ", " << ss.second << ") ";
    cout << endl;
    cout << "\tTarget: " << B.target << endl;
    cout << "\tIs leaf? " << B.isLeaf << endl;

    cout << endl;

    return;
}



// these are global vectors 
vector<bool> visited;
vector<int> parent;
vector<int> tree_stack;
vector<int> disc;
vector<int> low;
int visit_time;
vector<bool> is_art_point;
vector<int> times_art_point; // counts the multiplicity with which a node is art point; this is needed for source multiplicity

void findBCCs(int u, BCC &B, vector<BCC> &BCC_vector, vector<int> &source_neighbors)
{   
    // if(DEBUG) cout << "We are at node " << u << endl<<flush;
    tree_stack.push_back(u);
    // Count of children in DFS Tree (for root)
    int children = 0;
 
    // Mark the current node as visited
    visited[u] = true;
 
    // Initialize discovery time and lowpoint value
    disc[u] = low[u] = ++visit_time;
    

    // if(DEBUG) { cout << "Edges incident to " << u << " are: ";
    // for (auto inc : B.edges[u])
    // {
    //     cout << inc << " ";
    // }
    // cout << endl; }
    

    // Go through all neighbors of u
    for (auto v : B.edges[u]) {
        // if v is not visited yet, then make it a child of u in DFS tree and recur for it
        // cout << "Considering " << v << " as neighbor of " << u << endl; 
        if (!visited[v]) {
            // if(DEBUG) cout << "Considering unvisited neighbor " << v << " of " << u << endl;
            parent[v] = u;
            children++;
            findBCCs(v, B, BCC_vector, source_neighbors);

            // Check if the subtree rooted with v has
            // a connection to one of the ancestors of u
            low[u] = min(low[u], low[v]);

            // if u is not root and low value of one of its child is more than discovery value of u, 
            // OR if u is the root (returning from a child of the root identifies a BCC), then we close a BCC
            if ((parent[u] != -1 && low[v] >= disc[u]) || parent[u] == -1){ // here is where I close my articulation point
                BCC current_BCC;
                current_BCC.target = u; // u is the target of the current BCC
                current_BCC.nodes.push_back(u);
                is_art_point[u] = true;
                times_art_point[u]++;

                if(DEBUG) {cout << "Closing art point " << u << " because of " << v << endl;
                cout << "Current stack: ";
                for (auto ss : tree_stack)
                    cout << ss << " ";
                cout << endl;}
                
                if(DEBUG) {cout << "Current art points: ";
                for(int i = 0; i < is_art_point.size(); i++)
                {
                    if(is_art_point[i])
                        cout << i << " ";
                }
                cout << endl;}
                

                // for the current BCC we need to unstack until v
                int x = tree_stack.back();
                while(x != v){ // If we are in the "right" BCC add it to the vector   
                    current_BCC.nodes.push_back(x);
                    if(DEBUG) cout << "Adding node " << x << " to current BCC"<< endl;
                    
                    int source_multiplicity = 0;

                    // check if x is an art point; in this case it is a source with multiplicity equal to the number of children bccs
                    if(is_art_point[x])
                        source_multiplicity += times_art_point[x];
                    
                    // if furthermore x is a neighbor of the source, then it has a further 1 in its multiplicity
                    if(find(source_neighbors.begin(), source_neighbors.end(), x) != source_neighbors.end())
                        source_multiplicity++;

                    // at this point, if the multiplicity is greater than zero then we add x as a source 
                    if(source_multiplicity>0)
                        current_BCC.sources.push_back(make_pair(x, source_multiplicity));

                    // vector<pair<int,int>>::iterator find_source = find_if(current_BCC.sources.begin(), current_BCC.sources.end(), [&x](const pair<int, int>& p) { return p.first == x; });
                        
                    // if(find_source != current_BCC.sources.end()) { // if source is already present, increase its multiplicity
                    //     if(DEBUG) cout << "multiplicity of " << find_source->first<<  " before: " << find_source->second << endl; 
                    //     find_source->second++;
                    //     if(DEBUG) cout << "multiplicity after: " << find_source->second << endl; 
                    // }
                    // else // otherwise, create new source of multiplicity 1
                    //     current_BCC.sources.push_back(make_pair(x, 1));
                    

                    tree_stack.pop_back();
                    x = tree_stack.back();
                }
                tree_stack.pop_back(); 

                if(x!= v) // always finish at v
                    throw logic_error("Destacked more than what we should have");

                if(DEBUG) cout << "Adding node " << v << " to current BCC"<< endl;

                // TODO: FOR THE MULTIPLICITY, WE NEED TO COUNT HOW MANY BCCS CONTAIN THE NODE!!!
                // again, add v to sources if necessary
                // if(is_art_point[v] || find(source_neighbors.begin(), source_neighbors.end(), v) != source_neighbors.end()){
                //     if(DEBUG) cout << "Adding node " << v << " to sources"<< endl;
                //     vector<pair<int,int>>::iterator find_source = find_if(current_BCC.sources.begin(), current_BCC.sources.end(), [&v](const pair<int, int>& p) { return p.first == v; });
                    
                //     if(find_source != current_BCC.sources.end()) { // if source is already present, increase its multiplicity
                //         if(DEBUG) cout << "multiplicity of " << find_source->first<<  " before: " << find_source->second << endl; 
                //         find_source->second++;
                //         if(DEBUG) cout << "multiplicity after: " << find_source->second << endl; 
                //     }
                //     else // otherwise, create new source of multiplicity 1
                //         current_BCC.sources.push_back(make_pair(v, 1));
                // }
                
                int source_multiplicity = 0;

                // check if x is an art point; in this case it is a source with multiplicity equal to the number of children bccs
                if(is_art_point[v])
                    source_multiplicity += times_art_point[v];
                
                // if furthermore x is a neighbor of the source, then it has a further 1 in its multiplicity
                if(find(source_neighbors.begin(), source_neighbors.end(), v) != source_neighbors.end())
                    source_multiplicity++;

                // at this point, if the multiplicity is greater than zero then we add x as a source 
                if(source_multiplicity>0)
                    current_BCC.sources.push_back(make_pair(v, source_multiplicity));

                current_BCC.nodes.push_back(v); 

                current_BCC.isLeaf = false;
                current_BCC.myID = ++maxID;
                current_BCC.prodAncestors = -1;
                current_BCC.num_edges = 0;

                // for each node x, iterate through incident edges for it in B and add edges if the other endpoint is also in currentBCC
                for(auto x: current_BCC.nodes) {
                    current_BCC.edges.insert({x, {}}); // start with empty vector for all nodes
                    for (auto y: B.edges.at(x))
                    {
                        if(find(current_BCC.nodes.begin(), current_BCC.nodes.end(), y) != current_BCC.nodes.end()){
                            current_BCC.edges[x].push_back(y); // add y to edges' vector of x
                            current_BCC.num_edges++;
                        }
                    }
                }

                current_BCC.num_edges = current_BCC.num_edges/2;

                // add the BCC to the final vector
                BCC_vector.push_back(current_BCC);

            }
        }

        // Update low value of u for parent function calls.
        else if (v != parent[u])
            low[u] = min(low[u], disc[v]);
    }

    // root is art point only if it has more than one child
    if(parent[u] == -1 && children == 1)
        is_art_point[u] = false;    
    
    return;
}



// INPUT: a BCC B and one of its sources s
// OUTPUT: a vector of non trivial BCCs which form the block-cut tree of B after the removal of s
vector<BCC> expand(BCC &B, int s){
    // check if s is a source of to B
    if(find_if(B.sources.begin(), B.sources.end(), [&s](const pair<int, int>& p) { return p.first == s; }) == B.sources.end())
        throw invalid_argument("Source provided is not a source for the BCC");

    // remove s from nodes, incident edges of s from edges
    B.nodes.erase(find(B.nodes.begin(), B.nodes.end(), s));
    vector<int> source_neighbors = B.edges[s]; // save neighbors of s for later
    B.edges.erase(B.edges.find(s));
    for (auto neigh: source_neighbors) // also delete from neighbors' list
    {
        B.edges[neigh].erase(find(B.edges[neigh].begin(), B.edges[neigh].end(), s));
    }
    

    if(DEBUG) {
        cout << "Erased source; current graph is " << endl;
        for(auto x : B.nodes){
            cout << x << ": ";
            for(auto y : B.edges[x]){
                cout << y << " ";
            }
            cout << endl<<flush;
        }

        cout << "Neighbors of the source were: ";
        for(auto neigh: source_neighbors)
            cout << neigh << " ";
        cout <<endl;
    }


    // initialize vectors (maybe make them global)
    tree_stack.erase(tree_stack.begin(), tree_stack.end());
    visit_time = 0;


    // only initialize BCC nodes
    for(auto u : B.nodes){
        visited[u] = false;
        parent[u] = -2;
        is_art_point[u] = false;
        times_art_point[u]=0;
    }
    
    parent[B.target] = -1;


    // start visit from target of BCC (we first close BCCs that contain neighbors of s)
    vector<BCC> decomposed_BCC = {};
    findBCCs(B.target, B, decomposed_BCC, source_neighbors);


    // after call for art points, IF STACK IS NONEMPTY finish unstacking and form last new BCC
    // note: by non-empty we mean at least two elements, as if t is an articulation point then it must remain in the stack
    if(tree_stack.size() > 1){
        if(DEBUG){
            cout << "Stack is non-empty after function; it is: ";
            for (auto ss : tree_stack)
                cout << ss << " ";
            cout << endl;
        }

        BCC current_BCC;
        current_BCC.target = B.target; // B.target is the target of the current BCC
        // is_art_point[B.target] = true;

        // for the current BCC we need to unstack until the stack is empty
        // since stack is actually a vector, go through the whole vector
        for(auto stack_el : tree_stack){
            current_BCC.nodes.push_back(stack_el);

            // if(is_art_point[stack_el] || find(source_neighbors.begin(), source_neighbors.end(), stack_el) != source_neighbors.end()){
            //     vector<pair<int,int>>::iterator find_source = find_if(current_BCC.sources.begin(), current_BCC.sources.end(), [&stack_el](const pair<int, int>& p) { return p.first == stack_el; });
                
            //     if(find_source != current_BCC.sources.end()){ // if source is already present, increase its multiplicity
            //             if(DEBUG) cout << "multiplicity of " << find_source->first<<  " before: " << find_source->second << endl; 
            //             find_source->second++;
            //             if(DEBUG) cout << "multiplicity after: " << find_source->second << endl; 
            //         }
            //     else // otherwise, create new source of multiplicity 1
            //         current_BCC.sources.push_back(make_pair(stack_el, 1));
            // }

            int source_multiplicity = 0;

            // check if x is an art point; in this case it is a source with multiplicity equal to the number of children bccs
            if(is_art_point[stack_el])
                source_multiplicity += times_art_point[stack_el];
            
            // if furthermore x is a neighbor of the source, then it has a further 1 in its multiplicity
            if(find(source_neighbors.begin(), source_neighbors.end(), stack_el) != source_neighbors.end())
                source_multiplicity++;

            // at this point, if the multiplicity is greater than zero then we add x as a source 
            if(source_multiplicity>0)
                current_BCC.sources.push_back(make_pair(stack_el, source_multiplicity));
        }

        current_BCC.isLeaf = false;
        current_BCC.myID = ++maxID;
        current_BCC.prodAncestors = -1;
        current_BCC.num_edges=0;

        // fill edges right away: for each node x, iterate through incident edges for it in B and add edges if the other endpoint is also in currentBCC
        for(auto x: current_BCC.nodes) {
            current_BCC.edges.insert({x, {}}); // start with empty vector for all nodes
            for (auto y: B.edges[x])
            {
                if(find(current_BCC.nodes.begin(), current_BCC.nodes.end(), y) != current_BCC.nodes.end()){
                    current_BCC.edges[x].push_back(y); // add y to edges' vector of x
                    current_BCC.num_edges++;
                }
            }
        }

        current_BCC.num_edges = current_BCC.num_edges/2;

        // add the BCC to the final vector
        decomposed_BCC.push_back(current_BCC);
    }

    // need to set up personal bounds, isLeaf parameter by iterating through all new components
    for (int comp = 0; comp < decomposed_BCC.size(); comp++)
    {   
        decomposed_BCC[comp].personalBound = decomposed_BCC[comp].num_edges - decomposed_BCC[comp].nodes.size() +1; // personal bound given by cyclomatic lemma
        
        // it is a leaf if the sources vector has size one, and the only node belongs to the neighbors of s
        if(decomposed_BCC[comp].sources.size() == 1 && find(source_neighbors.begin(), source_neighbors.end(), decomposed_BCC[comp].sources[0].first)!= source_neighbors.end())
            decomposed_BCC[comp].isLeaf = true;

        // if(DEBUG) cout << "Set up for component number " << comp << endl;
        // if(DEBUG) cout << "Number of nodes: " << decomposed_BCC[comp].nodes.size()<<endl;
        // if(DEBUG) cout << "Number of edges: " << decomposed_BCC[comp].num_edges << endl;
        // if(DEBUG) cout << "Number of sources: " << decomposed_BCC[comp].sources.size()<<endl;
        // if(DEBUG) cout << "Personal bound is: " << decomposed_BCC[comp].personalBound << endl;
        // if(DEBUG) cout << "Component is leaf? " << decomposed_BCC[comp].isLeaf << endl;
    }
    return decomposed_BCC;
}


int main(int argc, char** argv) { 
    #ifdef RANDOM_NODE_EXTRACTION
        srand(time(NULL));
    #endif

    BCC G;
    N = 13;
    maxID = 1;

    visited.resize(N);
    parent.resize(N);
    disc.resize(N);
    low.resize(N);
    is_art_point.resize(N);
    times_art_point.resize(N);

    // G.sources = {0};
    // G.target = 10;
    // G.myID = 1;
    // G.nodes = {0,1,2,3,4,5,6,7,8,9,10};
    // G.edges.insert({0, {1,4,6,7}});
    // G.edges.insert({1, {0,2,3}});
    // G.edges.insert({2, {1,3}});
    // G.edges.insert({3, {1,2,4,5,6}});
    // G.edges.insert({4, {3,0,6}});
    // G.edges.insert({5, {3,6}});
    // G.edges.insert({6, {3,4,5,0,7,8}});
    // G.edges.insert({7, {6,8,0}});
    // G.edges.insert({8, {7,6,10,9}});
    // G.edges.insert({9, {8,10}});
    // G.edges.insert({10, {8,9}});

    G.sources = {make_pair(0, 1)};
    G.target = 12;
    G.myID = 1;
    G.nodes = {0,1,2,3,4,5,6,7,8,9,10,11,12};
    G.edges.insert({0, {2,7,8,10,11}});
    G.edges.insert({1, {2,3}});
    G.edges.insert({2, {1,0,3}});
    G.edges.insert({3, {1,2,4,5}});
    G.edges.insert({4, {3,5,12}});
    G.edges.insert({5, {4,3,6,7,12}});
    G.edges.insert({6, {5,7}});
    G.edges.insert({7, {5,6,0}});
    G.edges.insert({8, {0,9,10}});
    G.edges.insert({9, {8,10}});
    G.edges.insert({10, {0,8,9,11,12}});
    G.edges.insert({11, {12,10,0}});
    G.edges.insert({12, {4,5,10,11}});

    cout << "Original graph: " << endl;
    for(auto x : G.nodes){
        cout << x << ": ";
        for(auto y : G.edges[x]){
            cout << y << " ";
        }
        cout << endl;
    }


    cout << "Expanding with respect to source s= 0"<< endl;
    vector<BCC> result = expand(G, 0);

    cout << endl << "Found " << result.size() << " biconnected components:" << endl; 

    for (int i = 0; i < result.size(); i++)
    {
        cout << "B_"  << i << ": ";
        printBCC(result[i]);
    }

    return 0;
}