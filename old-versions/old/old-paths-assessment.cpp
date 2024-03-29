#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>
#include <chrono>
#include <cstdint>
#include <queue>

uint64_t timeMs() {
  using namespace std::chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

#define DEBUG true
// #define DEBUG false
#define RANDOM_NODE_EXTRACTION

using namespace std;

int N; // size of the graph/biggest component; used for global vectors
int maxID = 0; // maximum ID of a BCC
int s, t; // original source and target for st-paths problem
int z; // value to assess for
long MAX_TIME;
unsigned long running_bound = 0;

// ----------- global vectors for findBCCs
vector<bool> visited;
vector<int> parent;
vector<int> tree_stack;
vector<int> disc;
vector<int> low;
int visit_time;
vector<bool> is_art_point;
vector<int> times_art_point; // counts the multiplicity with which a node is art point; this is needed for source multiplicity
//-------------

struct BCC {
    vector<pair<int,int>> sources; // vector of paths' sources contained in the BCC. Multiplicity will be -1 for sources which are ONLY art points 
    int target; // target of the BCC = only articulation point leading to t 
    long personalBound; // personal bound of current BCC = |edges| - |nodes| + 1
    long prodAncestors; // product of the personal bounds of all BCCs in the bead string from current to root
    int myID; // personal ID of the BCC
    bool isLeaf; // true iff the BCC is a leaf of the block-cut tree

    vector<int> nodes; // vector of node indices forming the current BCC
    unordered_map<int, vector<int>> edges; // unordered map assigning to every node its incident edges
    int num_edges;
};


// node of the multi-source tree data structure
struct mptnode {
    int corrID; // ID of corresponding BCC
    BCC* corrBCC; // pointer to graph of corresponding BCC

    vector<mptnode*> children; // vector of pointers to its children
    mptnode* parent; // pointer to parent node
};


vector<mptnode*> tree_leaves; // vector of pointers to leaves of the multisource paths tree  
mptnode* mptree; // global structure; this is a pointer to the root of the multi-source paths tree


// vector<BCC*> tree_leaves; // vector of pointers to BCC leaves of the multisource paths tree  
// vector<BCC*> all_tree_nodes; // all the BCCs composing the multisource paths tree



// create graph from NDE file filename (WE REMOVE MULTI-EDGES AND SELF LOOPS)
void create_graph(char* filename)
{
    FILE* input_graph = fopen(filename, "r");

    // file contains number of nodes, number of edges at fist line
    // and then one edge per line
    fscanf(input_graph, "%d", &N);

    // if(DEBUG) cout << "Input graph has " << N << " nodes."<< endl;

    if(s >= N || t >= N){
        cout << "Input graph has " << N << " nodes. "<< endl;
        throw invalid_argument("Source or target are out of bounds");
    }

    visited.resize(N);
    parent.resize(N);
    disc.resize(N);
    low.resize(N);
    is_art_point.resize(N);
    times_art_point.resize(N);

    BCC* G = new BCC;
    G->sources = {make_pair(s, 1)};
    G->target = t;
    G->personalBound = 1;
    G->prodAncestors = 1;
    G->myID = maxID++;
    G->isLeaf = true;
    G->num_edges = 0;

    
    int u, v;
    // we need to skip the first N rows, plus add 0,...,N to nodes
    for(int i = 0; i < N; i++){
        fscanf(input_graph, "%d %d", &u, &v);
        G->nodes.push_back(i);
    }
    
    // for(int i=0; i<M; i++)
    while(fscanf(input_graph,"%d %d",&u, &v) != EOF)
    {
        if(DEBUG) cout << "Considering pair " << u << " " << v << endl;
        // make sure no self-loops or multiedges are created 
        if (u != v && find(G->edges[u].begin(), G->edges[u].end(), v) == G->edges[u].end()){
            G->edges[u].push_back(v);
            G->edges[v].push_back(u);
            G->num_edges++;
        }   
    }

    cout << "Input graph has " << N << " nodes and " << G->num_edges << " edges. "<< endl;

    fclose(input_graph);

    // G->personalBound = G->num_edges - G->nodes.size() + 1;// NO!!! ONLY TRUE IF BICONNECTED

    // tree_leaves.push_back(G); // initialize the leaf nodes and all nodes
    // all_tree_nodes.push_back(G);
    // running_bound = G->personalBound;

    // initialize multisource paths tree as root node with single child G
    mptnode* Gnode = new mptnode;
    Gnode->corrBCC = G;
    Gnode->corrID = G->myID;
    Gnode->parent = mptree;
    Gnode->children = {};

    mptree = new mptnode; 
    mptree->corrBCC = NULL;
    mptree->corrID = -1;
    mptree->parent = NULL;
    mptree->children = {Gnode};


    tree_leaves.push_back(Gnode); // initialize the good leaves as the root of the tree

    return;
}
 


void printBCC(BCC* B){
    cout << "\tGraph: ";
    for(auto x : B->nodes){
        cout <<"\t"<< x << ": [";
        for(auto y : B->edges[x]){
            cout << y << " ";
        }
        cout << "]"; 
    }
    cout << endl;
    cout << "\tPersonal bound: " << B->personalBound << endl;
    cout << "\tProduct of ancestors: " << B->prodAncestors << endl;
    cout << "\tID: " << B->myID << endl;
    cout << "\tSources: ";
    for(auto ss : B->sources)
        cout << "(" << ss.first << ", " << ss.second << ") ";
    cout << endl;
    cout << "\tTarget: " << B->target << endl;
    cout << "\tIs leaf? " << B->isLeaf << endl;

    cout << endl;

    return;
}

// print all my children, then endl and myself
void printMPtree(mptnode * tree_node){
    for (auto child : tree_node->children){
        printMPtree(child);
        cout << endl;    
    }
    cout << tree_node->corrID << " ";

    return;
}

// delete all data structures for the algorithm
void delete_all(mptnode* tree_node){
    for (auto child : tree_node->children)
        delete_all(child);
    
    delete tree_node->corrBCC;
    delete tree_node;

    return;
}

void findBCCs(int u, BCC* B, vector<BCC*> &BCC_vector, vector<int> &source_neighbors, int og_multiplicity)
{   
    // if(DEBUG) cout << "We are at node " << u << endl<<flush;
    tree_stack.push_back(u);
    int children = 0; // Count of children in DFS Tree (for root)
    visited[u] = true; // Mark the current node as visited
    disc[u] = low[u] = ++visit_time; // Initialize discovery time and lowpoint value
    
    // if(DEBUG) { cout << "Edges incident to " << u << " are: ";
    // for (auto inc : B->edges[u])
    // {
    //     cout << inc << " ";
    // }
    // cout << endl; }

    // Go through all neighbors of u
    for (auto v : B->edges[u]) {
        // if v is not visited yet, then make it a child of u in DFS tree and recur for it
        if (!visited[v]) {
            // if(DEBUG) cout << "Considering unvisited neighbor " << v << " of " << u << endl;
            parent[v] = u;
            children++;
            findBCCs(v, B, BCC_vector, source_neighbors, og_multiplicity);
            // if (DEBUG) cout << "Exited from recursive call for " << v << endl;
            // if (DEBUG) cout << "Values at this point: low[u]=" << low[u] << ", low[v]=" << low[v] << ", parent[u] =" << parent[u] << ", disc[u] =" << disc[u]<< endl;
            
            // Check if the subtree rooted with v has a connection to one of the ancestors of u
            low[u] = min(low[u], low[v]);

            // if u is not root and low value of one of its child is more than discovery value of u, 
            // OR if u is the root (returning from a child of the root identifies a BCC), then we close a BCC
            if ((parent[u] != -1 && low[v] >= disc[u]) || parent[u] == -1){ // here is where I close my articulation point
                if(DEBUG) {
                    cout << "Closing art point " << u << " because of " << v << endl;
                    cout << "Current stack: ";
                    for (auto ss : tree_stack)
                        cout << ss << " ";
                    cout << endl;
                }

                BCC* current_BCC = new BCC;
                current_BCC->target = u; // u is the target of the current BCC
                current_BCC->nodes = {};
                current_BCC->nodes.push_back(u);
                is_art_point[u] = true;
                times_art_point[u]++;
                
                if(DEBUG) {cout << "Current art points: "<<flush;
                for(int i = 0; i < is_art_point.size(); i++)
                {
                    if(is_art_point[i])
                        cout << i << " ";
                }
                cout << endl;}
                
                // for the current BCC we need to unstack until v
                int x = tree_stack.back();
                while(x != v){ // If we are in the "right" BCC add it to the vector   
                    current_BCC->nodes.push_back(x);
                    // if(DEBUG) cout << "Adding node x=" << x << " to current BCC"<< endl;

                    // if (DEBUG){
                    //     cout << "Stack is: ";
                    //     for (auto stack_el : tree_stack)
                    //         cout << stack_el << " ";
                    //     cout << endl;
                    // }
                    
                    
                    // multiplicity of the node is -1 if it is an art point; otherwise it is the multiplicity
                    // of the just-removed source summed to the possible multiplicity x could have, if it was a source of B
                    int source_multiplicity = 0;
                    bool art_pt = is_art_point[x];
                    vector<pair<int,int>>::iterator find_if_source = find_if(B->sources.begin(), B->sources.end(), [&x](const pair<int,int> &p){return p.first==x;});
                    bool old_source = (find_if_source != B->sources.end());
                    bool source_neigh = (find(source_neighbors.begin(), source_neighbors.end(), x) != source_neighbors.end());

                    if(old_source)
                        source_multiplicity+=find_if_source->second;
                    
                    if(source_neigh)
                        source_multiplicity+=og_multiplicity;
                    
                    if(source_multiplicity==0 && art_pt) // if none of the previous, but it is an art point, it still needs to be added as source with  mult -1
                        source_multiplicity = -1;
                        
                    // at this point, if the multiplicity is different from zero then we add x as a source 
                    if(source_multiplicity!=0)
                        current_BCC->sources.push_back(make_pair(x, source_multiplicity));

                    tree_stack.pop_back();
                    x = tree_stack.back();
                }
                tree_stack.pop_back(); 

                if(x!= v) // always finish at v
                    throw logic_error("Destacked more than what we should have");

                // if(DEBUG) cout << "Adding node v=" << v << " to current BCC"<< endl;
                
                // multiplicity of the node is -1 if it is an art point; otherwise it is the multiplicity
                // of the just-removed source summed to the possible multiplicity v could have, if it was a source of B
                int source_multiplicity = 0;
                bool art_pt = is_art_point[v];
                vector<pair<int,int>>::iterator find_if_source = find_if(B->sources.begin(), B->sources.end(), [&v](const pair<int,int> &p){return p.first==v;});
                bool old_source = (find_if_source != B->sources.end());
                bool source_neigh = (find(source_neighbors.begin(), source_neighbors.end(), v) != source_neighbors.end());

                if(old_source)
                    source_multiplicity+=find_if_source->second;
                
                if(source_neigh)
                    source_multiplicity+=og_multiplicity;
                
                if(source_multiplicity==0 && art_pt) // if none of the previous, but it is an art point, it still needs to be added as source with  mult -1
                    source_multiplicity = -1;
                    
                // at this point, if the multiplicity is different from zero then we add v as a source 
                if(source_multiplicity!=0)
                    current_BCC->sources.push_back(make_pair(v, source_multiplicity));


                current_BCC->nodes.push_back(v); 

                current_BCC->isLeaf = false;
                current_BCC->myID = ++maxID;
                current_BCC->prodAncestors = -1;
                current_BCC->num_edges = 0;

                // for each node x, iterate through incident edges for it in B and add edges if the other endpoint is also in currentBCC
                for(auto x: current_BCC->nodes) {
                    current_BCC->edges.insert({x, {}}); // start with empty vector for all nodes
                    for (auto y: B->edges.at(x))
                    {
                        if(find(current_BCC->nodes.begin(), current_BCC->nodes.end(), y) != current_BCC->nodes.end()){
                            current_BCC->edges[x].push_back(y); // add y to edges' vector of x
                            current_BCC->num_edges++;
                        }
                    }
                }
                current_BCC->num_edges = current_BCC->num_edges/2;
                // PERSONAL BOUND MUST BE CYCLOMATIC MULTIPLIED FOR THE NUMBER OF "VALID" SOURCES = NOT CONTAINING -1
                int cyclomatic_bound = current_BCC->num_edges - current_BCC->nodes.size() + 1; // initialize personal bound
                // current_BCC->personalBound = cyclomatic_bound;
                for (auto spair : current_BCC->sources)
                {
                    if(spair.second != -1)
                        current_BCC->personalBound += cyclomatic_bound*spair.second;
                }

                // also, set up if it is a leaf: that is, if it ONLY has sources with multiplicity > 0, that is, they are neighbors of the source
                if(find_if(current_BCC->sources.begin(), current_BCC->sources.end(), [](const pair<int, int>& p){return p.second == -1;}) == current_BCC->sources.end())
                    current_BCC->isLeaf = true;

                // add the BCC to the final vector
                BCC_vector.push_back(current_BCC);
            }
        } 
        else if (v != parent[u]) // Update low value of u for parent function calls.
            low[u] = min(low[u], disc[v]);
    }

    // root is art point only if it has more than one child
    if(parent[u] == -1 && children == 1)
        is_art_point[u] = false;    
    
    return;
}


// INPUT: BCC newB that has been exploded into decomposed_BCC vector of BCC, pointer to the root of the new tree for this vector of BCC and to node_to_replace, node of the tree to be replaced with the new subtree
// OUTPUT: Creates the full multisource paths tree for decomposed_BCC with root new_tree, also filling the corresponding prodAncestor bounds and adding leaves to leaf vector 
void update_tree(BCC* newB, vector<BCC*> &decomposed_BCC, mptnode* node_to_replace){
    // if(DEBUG) cout << "Updating BCC with ID " << newB->myID << "; it has target " << newB->target << endl;
    // if(DEBUG) cout << "Corresponding node to replace is the correct one? " << (node_to_replace->corrBCC == newB) << "; it has target " << node_to_replace->corrBCC->target << endl;
    int final_target = newB->target;
    vector<mptnode*> root_nodes; // roots of the new subtree
    vector<mptnode*> tree_nodes_stack;
    mptnode* new_subtree;

    //  PROBLEM: SEVERAL CAN HAVE TARGET EQUAL TO THE FINAL ONE
    vector<BCC*>::iterator find_root = find_if(decomposed_BCC.begin(), decomposed_BCC.end(), [&final_target](const BCC* comp){return comp->target == final_target;});
    if(find_root == decomposed_BCC.end())
        throw logic_error("Error when exploding: cannot find component containing target");

    while (find_root!=decomposed_BCC.end())
    {
        find_root = find_if(find_root, decomposed_BCC.end(), [&final_target](const BCC* comp){return comp->target == final_target;});
        if(find_root!= decomposed_BCC.end()){ // in this case, we found a component which has target equal to the current source of current_node
            // create a new node attached to current_node, push it to the tree_nodes_stack, and update the corresponding prodAncestors bound accordingly
            new_subtree =  new mptnode;
            new_subtree->corrBCC = *find_root;
            new_subtree->corrID = (*find_root)->myID;
            new_subtree->children = {};
            new_subtree->parent = node_to_replace->parent; // attach new tree to the parent of the leaf we exploded
            
            if(new_subtree->parent == mptree && node_to_replace->corrBCC->target != t)
                throw logic_error("Root of the multisource paths tree does not contain target");

            new_subtree->corrBCC->prodAncestors = node_to_replace->corrBCC->prodAncestors; // the product of the ancestors of the root is equal to the one of the old leaf node
            if(DEBUG) cout << "Product of ancestors of root node is " << new_subtree->corrBCC->prodAncestors<<endl;

            tree_nodes_stack.push_back(new_subtree);
            root_nodes.push_back(new_subtree);
            if(DEBUG) cout << "Adding component with ID="<< new_subtree->corrID << " to the tree node stack and to tree roots"<< endl;

            find_root++;
        }
    }  
    
    if(DEBUG){
        cout << "BCC IDs are: ";
        for (auto cc : decomposed_BCC){
            cout << cc->myID << " ";
        }
        cout << endl;
    }

    
    mptnode *current_node, *new_node;
    while (!tree_nodes_stack.empty()){
        current_node = tree_nodes_stack.back();
        tree_nodes_stack.pop_back();
        if(!current_node->corrBCC->isLeaf){ // only explore children if the corresponding BCC is not a leaf
            if (DEBUG) cout << "Considering non-leaf BCC with ID=" << current_node->corrBCC->myID << endl;
            for (auto spair : current_node->corrBCC->sources){ // go through all sources of the component 
                if (DEBUG) cout << "Considering source " << spair.first << endl;
                auto i = decomposed_BCC.begin();
                auto end = decomposed_BCC.end();
                while (i!=end)
                {
                    i = find_if(i, end, [&](const BCC* comp){return comp->target == spair.first;});
                    if(i!= end){ // in this case, we found a component which has target equal to the current source of current_node
                        // create a new node attached to current_node, push it to the tree_nodes_stack, and update the corresponding prodAncestors bound accordingly
                        if (DEBUG) cout << "Found BCC " << decomposed_BCC[i-decomposed_BCC.begin()]->myID << " to have target equal to "<< spair.first << endl;
                        new_node = new mptnode;
                        new_node->corrBCC = *i;
                        new_node->corrID = new_node->corrBCC->myID;
                        new_node->parent = current_node;
                        current_node->children.push_back(new_node);
                        new_node->corrBCC->prodAncestors = current_node->corrBCC->prodAncestors*current_node->corrBCC->personalBound;

                        tree_nodes_stack.push_back(new_node);
                        if(DEBUG) cout << "Adding component with ID="<< new_node->corrID << " to the tree node stack "<< endl;
                        i++;
                    }
                }  
            }
        }
        else{ // if current node is indeed a leaf, add it to the leaf vector of tree nodes
            if (DEBUG) cout << "Found leaf node for component with ID=" << current_node->corrBCC->myID << endl;
            tree_leaves.push_back(current_node);
        }
    }

    // if(DEBUG) cout << "Finished filling new tree" << endl;

    // any new root node is already attached to node_to_replace; we just need to remove the latter altogether
    if(root_nodes[0]->parent != NULL){
        if(DEBUG) cout << "New root has non-NULL parent" << endl;
        auto child_iterator = find(node_to_replace->parent->children.begin(), node_to_replace->parent->children.end(), node_to_replace);
        if(child_iterator == node_to_replace->parent->children.end())
            throw logic_error("Parent of node to replace does not have him as one of its children");
        node_to_replace->parent->children.erase(child_iterator); // remove the node as a child of its parent
    }
    else{ // if we were at the root, we need to assign mptree to the root of the new subtree
        mptree = new_subtree;
    }

    if(node_to_replace->children.size()!= 0) // extra check
        throw logic_error("Node to replace has children");

    // remove the full node and its BCC
    delete node_to_replace->corrBCC;
    delete node_to_replace;

    return;
}

// INPUT: a pointer to the leaf BCC we wish to expand
// OUTPUT: a vector of non trivial BCCs which form the block-cut tree of B after the removal of s
// void explode(BCC* B){
void explode(mptnode* leaf_node){
    BCC* B = leaf_node->corrBCC;

    // check if the node is indeed a leaf
    if(!B->isLeaf || leaf_node->children.size() > 0)
        throw invalid_argument("Trying to expand a non-leaf node of the multi-source paths tree");

    // choose a source at random 
    int sIndex = rand() % B->sources.size();
    pair<int,int> sBpair = B->sources[sIndex];

    // remember multiplicity of s, as it will be the starting one of all its neighbors
    int og_multiplicity = sBpair.second;
    int sB = sBpair.first;

    if(og_multiplicity == -1) // the BCC had children BCCs
        throw logic_error("Trying to expand a non-leaf node of the multi-source paths tree");


    // we need to consider two cases according to the size of sBpair:
    // 1) if it is 1, we remove B and substitute it with the new tree of components after the removal of s
    // 2) if it is >1, B stays and we modify a copy newB, which will be a new child of the parent of leaf_node
    // the same for the tree node
    // in any case, the selected source must be removed from the graph's sources
    if(sBpair.second > 1)
        B->sources[sIndex].second--;
    else
        B->sources.erase(B->sources.begin() + sIndex);

    BCC* newB = B; // start with a pointer to B
    mptnode* node_to_replace = leaf_node;
    if(B->sources.size()>0){ // if the number of sources is greater than zero (sB has been already removed), create an actual copy
        newB = new BCC;
        newB->sources = B->sources;
        newB->target = B->target;
        newB->personalBound = B->personalBound;
        newB->prodAncestors = B->prodAncestors;
        newB->myID = B->myID;
        newB->isLeaf = B->isLeaf;
        newB->nodes =B->nodes;
        newB->edges = B->edges;
        newB->num_edges = B->num_edges;

        node_to_replace = new mptnode;
        node_to_replace->corrID = leaf_node->corrID;
        node_to_replace->corrBCC = newB;
        node_to_replace->children = {};
        node_to_replace->parent = leaf_node->parent;
    }
    // now newB, node_to_replace hold the correct graph to modify

    // remove s from nodes, incident edges of s from edges
    newB->nodes.erase(find(newB->nodes.begin(), newB->nodes.end(), sB));
    vector<int> source_neighbors = newB->edges[sB]; // save neighbors of s for later

    newB->edges.erase(newB->edges.find(sB));
    for (auto neigh: source_neighbors) // also delete from neighbors' list
        newB->edges[neigh].erase(find(newB->edges[neigh].begin(), newB->edges[neigh].end(), sB));

    if(DEBUG) {
        cout << "Erased source; current graph is " << endl;
        for(auto x : newB->nodes){
            cout << x << ": ";
            for(auto y : newB->edges[x]){
                cout << y << " ";
            }
            cout << endl<<flush;
        }

        cout << "Neighbors of the source were: ";
        for(auto neigh: source_neighbors)
            cout << neigh << " ";
        cout <<endl;
    }

    // initialize support vectors for findBCC
    tree_stack.erase(tree_stack.begin(), tree_stack.end());
    visit_time = 0;
    // only initialize BCC nodes
    for(auto u : newB->nodes){
        visited[u] = false;
        parent[u] = -2;
        is_art_point[u] = false;
        times_art_point[u]=0;
    }
    parent[newB->target] = -1;

    // start visit from target of BCC (we first close BCCs that contain neighbors of s)
    vector<BCC*> decomposed_BCC = {}; // DO WE ALSO NEED TO BUILD A VECTOR OF TREE NODES???
    findBCCs(newB->target, newB, decomposed_BCC, source_neighbors, og_multiplicity);
    int og_leaves_size = tree_leaves.size(); // recall the size before adding new leaves


    if(DEBUG){
        cout << "BCC found are: "<<  endl;
        for (int i = 0; i < decomposed_BCC.size(); i++){
            cout << "B_"  << i << ": ";
            printBCC(decomposed_BCC[i]);
        }
    }

    // after call for art points, IF STACK IS NONEMPTY finish unstacking and form last new BCC
    // note: by non-empty we mean at least two elements, as if t is an articulation point then it must remain in the stack
    if(tree_stack.size() > 1){
        if(DEBUG){
            cout << "Stack is non-empty after function; it is: ";
            for (auto ss : tree_stack)
                cout << ss << " ";
            cout << endl;
        }

        BCC* current_BCC = new BCC;
        current_BCC->target = newB->target; // B.target is the target of the current BCC

        // for the current BCC we need to unstack until the stack is empty
        // since stack is actually a vector, go through the whole vector
        for(auto stack_el : tree_stack){
            current_BCC->nodes.push_back(stack_el);
            // multiplicity of the node is -1 if it is an art point; otherwise it is the multiplicity
            // of the just-removed source summed to the possible multiplicity x could have, if it was a source of B
            int source_multiplicity = 0;
            bool art_pt = is_art_point[stack_el];
            vector<pair<int,int>>::iterator find_if_source = find_if(B->sources.begin(), B->sources.end(), [&stack_el](const pair<int,int> &p){return p.first==stack_el;});
            bool old_source = (find_if_source != B->sources.end());
            bool source_neigh = (find(source_neighbors.begin(), source_neighbors.end(), stack_el) != source_neighbors.end());

            if(old_source)
                source_multiplicity+=find_if_source->second;
            
            if(source_neigh)
                source_multiplicity+=og_multiplicity;
            
            if(source_multiplicity==0 && art_pt) // if none of the previous, but it is an art point, it still needs to be added as source with  mult -1
                source_multiplicity = -1;
                
            // at this point, if the multiplicity is different from zero then we add x as a source 
            if(source_multiplicity!=0)
                current_BCC->sources.push_back(make_pair(stack_el, source_multiplicity));

        }

        current_BCC->isLeaf = false;
        current_BCC->myID = ++maxID;
        current_BCC->prodAncestors = -1;
        current_BCC->num_edges=0;

        // fill edges right away: for each node x, iterate through incident edges for it in B and add edges if the other endpoint is also in currentBCC
        for(auto x: current_BCC->nodes) {
            current_BCC->edges.insert({x, {}}); // start with empty vector for all nodes
            for (auto y: newB->edges[x])
            {
                if(find(current_BCC->nodes.begin(), current_BCC->nodes.end(), y) != current_BCC->nodes.end()){
                    current_BCC->edges[x].push_back(y); // add y to edges' vector of x
                    current_BCC->num_edges++;
                }
            }
        }

        current_BCC->num_edges = current_BCC->num_edges/2;

        // PERSONAL BOUND MUST BE CYCLOMATIC MULTIPLIED FOR THE NUMBER OF "VALID" SOURCES = NOT CONTAINING -1
        int cyclomatic_bound = current_BCC->num_edges - current_BCC->nodes.size() + 1; // initialize personal bound
        // current_BCC->personalBound = cyclomatic_bound;
        for (auto spair : current_BCC->sources)
        {
            if(spair.second != -1)
                current_BCC->personalBound += cyclomatic_bound*spair.second;
        }
        

        // also, set up if it is a leaf: that is, if it ONLY has sources with multiplicity > 0, that is, they are neighbors of the source
        if(find_if(current_BCC->sources.begin(), current_BCC->sources.end(), [](const pair<int, int>& p){return p.second == -1;}) == current_BCC->sources.end()){
            current_BCC->isLeaf = true;
        }

        // add the BCC to the final vector
        decomposed_BCC.push_back(current_BCC);
    }

    // need to build tree and update prodAncestors fields 
    mptnode* new_subtree; // this is a pointer to the root of the new tree; the root corresponds to the only component holding newB->target as target
    update_tree(newB, decomposed_BCC, node_to_replace); // MISSING: UPDATE RUNNING BOUND!!!!!

    // if(DEBUG) cout << "Exited tree update procedure" <<endl<<flush;
    if(DEBUG){
        cout << "After updating the tree the new BCCs are: "<<  endl;
        for (int i = 0; i < decomposed_BCC.size(); i++){
            cout << "B_"  << i << ": ";
            printBCC(decomposed_BCC[i]);
        }
    }

    return;
}



int main(int argc, char** argv) { 
    if(argc < 5){
        cout << "USAGE: " << argv[0] << " <graph-filename> <source> <target> <z> [MAX_TIME]\n";
        return 0;
    }


    if(argc >= 6) MAX_TIME = atoi(argv[5]);

    s = atoi(argv[2]);
    t = atoi(argv[3]);
    z = atoi(argv[4]);

    if(s==t)
        throw invalid_argument("Source and target must be different");
    if(z<0)
        throw invalid_argument("Parameter z must be a positive integer");

    create_graph(argv[1]);
    
    #ifdef RANDOM_NODE_EXTRACTION
        srand(time(NULL));
    #endif


    if(DEBUG){
        cout << "Original graph: " << endl<<flush;
        for(auto x : mptree->corrBCC->nodes){
            cout << x << ": ";
            for(auto y : mptree->corrBCC->edges[x]){
                cout << y << " ";
            }
            cout << endl;
        }
    }


    // cout << "Expanding with respect to source s="<< s << endl;
    explode(mptree);
    
    printMPtree(mptree);
    cout << endl;

    delete_all(mptree);

    return 0;
}