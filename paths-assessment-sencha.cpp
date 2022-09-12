#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>
#include <chrono>
#include <cstdint>
#include <queue>
#include <unordered_set>

uint64_t timeMs() {
  using namespace std::chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

// #define DEBUG true
#define DEBUG false


using namespace std;

int N; // size of the graph/biggest component; used for global vectors
int M; // number of edges of the graph
int maxID = 0; // maximum ID of a BCC
int s, t; // original source and target for st-paths problem
int z; // value to assess for
long MAX_TIME = -1;
unsigned long running_bound = 1;
unsigned long calls_performed = 0;

// ----------- global vectors for findBCCs
vector<bool> visited;
vector<int> parent;
vector<int> tree_stack;
vector<int> disc;
vector<int> low;
int visit_time;
vector<bool> is_art_point;
vector<int> times_art_point; // counts the multiplicity with which a node is art point; this is needed for source multiplicity
// -----------

struct BCC {
    vector<pair<int,int>> sources; // vector of paths' sources contained in the BCC. Multiplicity will be -1 for sources which are ONLY art points 
    int target; // target of the BCC = only articulation point leading to t 
    long personalBound; // personal bound of current BCC = |edges| - |nodes| + 1
    int bound_multiplier; // multiplier given by the number of sources and their multiplicity; this influences the whole leaf-to-root path
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
vector<mptnode*> tree_semi_leaves; // vector of pointers to nodes of mptree whose corresponding BCC is a leaf (i.e. at least one non-art point leaf), needed for bound update
mptnode* mptree; // global structure; this is a pointer to the root of the multi-source paths tree


// vector<BCC*> tree_leaves; // vector of pointers to BCC leaves of the multisource paths tree  
// vector<BCC*> all_tree_nodes; // all the BCCs composing the multisource paths tree

long inline leaf_to_root_bound(mptnode* tree_node){
    return tree_node->corrBCC->personalBound*tree_node->corrBCC->prodAncestors*tree_node->corrBCC->bound_multiplier;
}

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
    G->bound_multiplier = 1;
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
        // if(DEBUG) cout << "Considering pair " << u << " " << v << endl;
        // make sure no self-loops or multiedges are created 
        if (u != v && find(G->edges[u].begin(), G->edges[u].end(), v) == G->edges[u].end()){
            G->edges[u].push_back(v);
            G->edges[v].push_back(u);
            G->num_edges++;
        }   
    }

    // cout << "Input graph has " << N << " nodes and " << G->num_edges << " edges. "<< endl;

    fclose(input_graph);

    // G->personalBound = G->num_edges - G->nodes.size() + 1;// NO!!! ONLY TRUE IF BICONNECTED

    // tree_leaves.push_back(G); // initialize the leaf nodes and all nodes
    // all_tree_nodes.push_back(G);
    // running_bound = G->personalBound;

    M= G->num_edges;


    // ADD NEW DUMMY NODE AND BUILD THE TREE
    BCC* dummyBCC = new BCC;
    dummyBCC->sources = {make_pair(t, -1)};
    dummyBCC->target = N;
    dummyBCC->personalBound = 1;
    dummyBCC->bound_multiplier = 0;
    dummyBCC->prodAncestors = 1;
    dummyBCC->myID = -1;
    dummyBCC->nodes = {t, N};
    dummyBCC->edges[t].push_back(N);
    dummyBCC->edges[N].push_back(t);    
    dummyBCC->num_edges = 1;
    
    t = N;

    // initialize multisource paths tree as dummy root with single child G 
    mptree = new mptnode; 
    mptree->corrBCC = dummyBCC;
    mptree->corrID = dummyBCC->myID;
    mptree->parent = NULL;
    mptree->children = {};

    mptnode* Gnode = new mptnode;
    Gnode->corrBCC = G;
    Gnode->corrID = G->myID;
    Gnode->parent = mptree;
    Gnode->children = {};

    mptree->children.push_back(Gnode);

    // do not initialize tree_leaves as G will be immediately exploded
    // tree_leaves.push_back(Gnode); // initialize the good leaves as the root of the tree

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
    cout << "\tBound multiplier: " << B->bound_multiplier << endl;
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
void printTreeNodes(mptnode * tree_node){
    cout << tree_node->corrID << "(" << tree_node->corrBCC->nodes.size() << " nodes) -> "; 
    for (auto child : tree_node->children){
        cout << child->corrID << " ";
    }
    cout << endl;

    for (auto child : tree_node->children){
        printTreeNodes(child);    
    }

    return;
}

void printTreeComponents(mptnode * tree_node){
    printBCC(tree_node->corrBCC);

    for (auto child : tree_node->children){
        printTreeComponents(child);    
    }

    return;
}

void printMPtree(mptnode * tree_node){
    printTreeNodes(tree_node);
    // cout << "Where the components are: "<< endl;
    // printTreeComponents(tree_node);    

    return;
}

void printLeaves(){
    cout << "IDs of current leaves: ";
    for (auto tnode : tree_leaves)
        cout << tnode->corrID << " ";
    
    cout << endl;
}

void printSemiLeaves(){
    cout << "IDs of current semi-leaves: ";
    for (auto tnode : tree_semi_leaves)
        cout << tnode->corrID << " ";
    
    cout << endl;
}

// delete a single mptnode
void inline delete_mptnode(mptnode* tree_node){
    delete tree_node->corrBCC;
    delete tree_node;
}

// delete all data structures for the algorithm
void delete_all(mptnode* tree_node){
    for (auto child : tree_node->children)
        delete_all(child);
    
    delete tree_node->corrBCC;
    delete tree_node;

    return;
}

// ----------------------------------------------------------------------------------------------------------


inline mptnode* pick_random_leaf(){
    return tree_leaves[rand() % tree_leaves.size()];
}


inline mptnode* pick_smallest_leaf(){
    int min_size = N+1;
    mptnode* smallest = tree_leaves[0];

    for(auto x : tree_leaves){
        if(x->corrBCC->nodes.size() < min_size){
            min_size = x->corrBCC->nodes.size();
            smallest = x;
            if(min_size==2)
                return smallest;
        }
    }

    return smallest;
}

inline mptnode* pick_biggest_leaf(){
    int max_size = 0;
    mptnode* biggest = tree_leaves[0];

    for(auto x : tree_leaves){
        if(x->corrBCC->nodes.size() > max_size){
            max_size = x->corrBCC->nodes.size();
            biggest = x;
        }
    }

    return biggest;
}


inline mptnode* pick_maxbound_leaf(){
    int max_bound = 0;
    mptnode* boundest = tree_leaves[0];

    for(auto x : tree_leaves){
        if(x->corrBCC->personalBound > max_bound){
            max_bound = x->corrBCC->personalBound;
            boundest = x;
        }
    }

    return boundest;
}


inline mptnode* pick_LIFO_leaf(){
    return tree_leaves[tree_leaves.size()-1];
}



void findBCCs(int u, BCC* B, vector<BCC*> &BCC_vector, vector<int> &source_neighbors, int og_multiplicity)
{   
    // if(DEBUG) cout << "We are at node " << u << endl<<flush;
    tree_stack.push_back(u);
    int children = 0; // Count of children in DFS Tree (for root)
    visited[u] = true; // Mark the current node as visited
    disc[u] = ++visit_time; // Initialize discovery time and lowpoint value
    low[u] = disc[u];
    
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
                // if(DEBUG) {
                //     cout << "Closing art point " << u << " because of " << v << endl;
                //     cout << "Current stack: ";
                //     for (auto ss : tree_stack)
                //         cout << ss << " ";
                //     cout << endl;
                // }

                unordered_set<int> curr_nodes = {};

                BCC* current_BCC = new BCC;
                current_BCC->target = u; // u is the target of the current BCC
                current_BCC->nodes = {};
                current_BCC->nodes.push_back(u);
                curr_nodes.insert(u);

                is_art_point[u] = true;
                times_art_point[u]++;
                
                // if(DEBUG) {cout << "Current art points: "<<flush;
                // for(int i = 0; i < is_art_point.size(); i++)
                // {
                //     if(is_art_point[i])
                //         cout << i << " ";
                // }
                // cout << endl;}
                
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
                    curr_nodes.insert(x);
                    
                    
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
                curr_nodes.insert(v);

                current_BCC->isLeaf = false;
                current_BCC->myID = ++maxID;
                current_BCC->prodAncestors = -1;
                current_BCC->num_edges = 0;
                current_BCC->bound_multiplier = 0;

                // for each node x, iterate through incident edges for it in B and add edges if the other endpoint is also in currentBCC
                for(auto x: current_BCC->nodes) {
                    
                    current_BCC->edges.insert({x, {}}); // start with empty vector for all nodes
                    for (auto y: B->edges[x])
                    {
                        // if(find(current_BCC->nodes.begin(), current_BCC->nodes.end(), y) != current_BCC->nodes.end()){
                        //     current_BCC->edges[x].push_back(y); // add y to edges' vector of x
                        //     current_BCC->num_edges++;
                        // }

                        if(curr_nodes.find(y)!= curr_nodes.end()){
                            current_BCC->edges[x].push_back(y); // add y to edges' vector of x
                            current_BCC->num_edges++;
                        }
                    }
                }
                current_BCC->num_edges = current_BCC->num_edges/2;
                // PERSONAL BOUND MUST BE CYCLOMATIC MULTIPLIED FOR THE NUMBER OF "VALID" SOURCES = NOT CONTAINING -1
                // int cyclomatic_bound = current_BCC->num_edges - current_BCC->nodes.size() + 1; // initialize personal bound
                current_BCC->personalBound = current_BCC->num_edges - current_BCC->nodes.size() + 1; 

                if(current_BCC->nodes.size() == 2) current_BCC->personalBound = 1;

                bool only_art_pts = true;
                
                // if(DEBUG) cout << "Number of nodes in this BCC is " << current_BCC->nodes.size()<<endl;
                // current_BCC->personalBound = cyclomatic_bound;
                for (auto spair : current_BCC->sources)
                {
                    if(spair.second != -1){
                        only_art_pts = false;
                        // current_BCC->personalBound += cyclomatic_bound*spair.second;
                        current_BCC->bound_multiplier += spair.second;
                    }
                }

                // if(only_art_pts) current_BCC->personalBound = cyclomatic_bound;

                // if(DEBUG) cout << "Personal bound is " << current_BCC->personalBound <<endl;

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
    // int og_leaves = tree_leaves.size(); // save for later vector<mptnode*> root_nodes; // roots of the new subtree
    // int og_semi_leaves = tree_semi_leaves.size();
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
            new_subtree->parent->children.push_back(new_subtree);

            // if(DEBUG) cout << "Just attached root of new tree with ID " << new_subtree->corrID << " to parent with ID " << new_subtree->parent->corrID << endl;
            
            if(new_subtree->parent == mptree && node_to_replace->parent->corrBCC->target != t)
                throw logic_error("Root of the multisource paths tree does not contain target");

            new_subtree->corrBCC->prodAncestors = node_to_replace->corrBCC->prodAncestors; // the product of the ancestors of the root is equal to the one of the old leaf node
            // if(DEBUG) cout << "Product of ancestors of root node is " << new_subtree->corrBCC->prodAncestors<<endl;

            tree_nodes_stack.push_back(new_subtree);
            root_nodes.push_back(new_subtree);
            // if(DEBUG) cout << "Adding component with ID="<< new_subtree->corrID << " to the tree node stack and to tree roots"<< endl;

            find_root++;
        }
    }  

    // vector<mptnode*> tree_nodes_stack = {new_subtree};
    mptnode *current_node, *new_node;
    while (!tree_nodes_stack.empty()){
        current_node = tree_nodes_stack.back();
        tree_nodes_stack.pop_back();
        // if (DEBUG)
        // {
        //     cout << "Considering component with ID " << current_node->corrID << endl;
        //     cout << "Sources are: ";
        //     for (auto spair : current_node->corrBCC->sources)
        //     {
        //         cout << "(" << spair.first << ", " <<spair.second << ")\t";
        //     }
        //     cout << endl;
        // }
        
        // if(!current_node->corrBCC->isLeaf){ // only explore children if the corresponding BCC is not a leaf
        // if(find_if(current_node->corrBCC->sources.begin(), current_node->corrBCC->sources.end(), [&](const pair<int,int> p){return p.second == -1;}) == current_node->corrBCC->sources.end()){
        //     if (DEBUG) cout << "Considering non-leaf BCC with ID=" << current_node->corrBCC->myID << endl;
        int children_count = 0; // keep count of how many real sources
        for (auto spair : current_node->corrBCC->sources){ // go through all sources of the component
            // if (DEBUG) cout << "Considering source " << spair.first << endl;
            auto i = decomposed_BCC.begin();
            auto end = decomposed_BCC.end();
            while (i!=end)
            {
                i = find_if(i, end, [&](const BCC* comp){return comp->target == spair.first;}); // this handles possibility of not having any child with that target
                if(i!= end){ // in this case, we found a component which has target equal to the current source of current_node
                    // create a new node attached to current_node, push it to the tree_nodes_stack, and update the corresponding prodAncestors bound accordingly
                    children_count++; // found a child
                    // if (DEBUG) cout << "Found BCC " << decomposed_BCC[i-decomposed_BCC.begin()]->myID << " to have target equal to "<< spair.first << endl;
                    new_node = new mptnode;
                    new_node->corrBCC = *i;
                    new_node->corrID = new_node->corrBCC->myID;
                    new_node->parent = current_node;
                    current_node->children.push_back(new_node);
                    // if(DEBUG){
                    //     cout << "current node prod ancestors: " << current_node->corrBCC->prodAncestors<<endl;
                    //     cout << "current node personal bound: " << current_node->corrBCC->personalBound<<endl;
                    // }
                    new_node->corrBCC->prodAncestors = current_node->corrBCC->prodAncestors*current_node->corrBCC->personalBound;

                    tree_nodes_stack.push_back(new_node);
                    // if(DEBUG) cout << "Adding component with ID="<< new_node->corrID << " to the tree node stack "<< endl;
                    i++;
                }
            }  
        }
        
        if(children_count==0){ // if current node is indeed a leaf, add it to the leaf vector of tree nodes
            // if (DEBUG) cout << "Found leaf node for component with ID=" << current_node->corrBCC->myID << endl;
            tree_leaves.push_back(current_node);
            // if(DEBUG) cout << "Size of tree leaves: " << tree_leaves.size() << endl;
        }

        if(children_count > 0 && current_node->corrBCC->isLeaf) { // if the corresponding BCC is leaf but it is not an actual tree leaf, then current node is a semi-leaf of the tree
            tree_semi_leaves.push_back(current_node);
        }
    }

    // if(DEBUG) cout << "Finished filling new tree" << endl;

    // new subtree is already attached to node_to_replace; we just need to remove the latter altogether
    if(new_subtree->parent != NULL){
        // if(DEBUG) cout << "new subtree has non-NULL parent" << endl;
        auto child_iterator = find(node_to_replace->parent->children.begin(), node_to_replace->parent->children.end(), node_to_replace);
        if(child_iterator == node_to_replace->parent->children.end())
            throw logic_error("Parent of node to replace does not have him as one of its children");
        node_to_replace->parent->children.erase(child_iterator); // remove the node as a child of its parent
    }
    else{ // if we were at the root, we need to assign mptree to the root of the new subtree
        // mptree = new_subtree;
    }

    if(node_to_replace->children.size()!= 0) // extra check
        throw logic_error("Node to replace has children");

    if(DEBUG){
        printLeaves();
        printSemiLeaves();
    }

    return;
}

// INPUT: a pointer to the leaf BCC we wish to expand
// OUTPUT: a vector of non trivial BCCs which form the block-cut tree of B after the removal of s
// void explode(BCC* B){
void explode(mptnode* leaf_node){
    BCC* B = leaf_node->corrBCC;

    // cout << "Considering node with ID " << B->myID << flush << endl;
    // printMPtree(mptree);
    // cout << "Finished print " << endl << flush;

    // treat case of trivial component (two nodes)
    if(B->nodes.size() == 2)
    { 
        if(DEBUG) cout << "Component is trivial! " << endl<< flush;

        // cout << "Component is trivial! " << endl<< flush;
        // printBCC(B);

        if(leaf_node->corrBCC->bound_multiplier != leaf_node->corrBCC->sources[0].second)
            throw logic_error("The bound multiplier of a trivial BCC is different from its source multiplicity");

        if(leaf_node->corrBCC->sources.size() > 1)
            throw logic_error("A trivial BCC has more than one source: impossible!");

        // running bound is unchanged as component personal bound was 1 

        // find the source of the parent that corresponds to the target of B
        mptnode * parent_node = leaf_node->parent;
        vector<pair<int,int>>::iterator find_source = find_if(parent_node->corrBCC->sources.begin(), parent_node->corrBCC->sources.end(), [&B](const pair<int,int> p){return p.first == B->target;});
        if(find_source == parent_node->corrBCC->sources.end())
            throw logic_error("Target of child node is not a source of parent node");

        // assign multiplicity of source of B to the find_source source and also add it to the bound multiplier of the parent component
        int sB_mult = B->sources[0].second;
        if(parent_node->corrBCC->sources[find_source-parent_node->corrBCC->sources.begin()].second == -1) // if it was art point, set its multiplicity to sB_mult
            parent_node->corrBCC->sources[find_source-parent_node->corrBCC->sources.begin()].second = sB_mult;
        else // otherwise, add the multiplicity to the one it already has
            parent_node->corrBCC->sources[find_source-parent_node->corrBCC->sources.begin()].second += sB_mult;

        // in any case, add it to the bound multiplier
        parent_node->corrBCC->bound_multiplier += sB_mult;

        // remove the leaf, and add the parent to leaves or semi-leaves according to whether it has more children
        vector<mptnode*>::iterator find_child = find(parent_node->children.begin(), parent_node->children.end(), leaf_node);
        if(find_child == parent_node->children.end())
            throw logic_error("Child is not in children vector of parent");
        
        parent_node->children.erase(find_child); // remove leaf node from parent's children

        // check if parent was a semi leaf
        vector<mptnode*>::iterator is_semi_leaf = find(tree_semi_leaves.begin(), tree_semi_leaves.end(), parent_node);

        if(parent_node->children.size() == 0){ // if leaf_node was the only child, then parent node becomes leaf node
            tree_leaves.push_back(parent_node);
            parent_node->corrBCC->isLeaf = true; // set the corresponding BCC as a leaf
            if(is_semi_leaf != tree_semi_leaves.end()) // if it was also a semi-leaf, it must be removed from that vector
                tree_semi_leaves.erase(is_semi_leaf);
        }       
        else{ // here it is surely a semi-leaf, so we add it to them if not already in the vector
            if(is_semi_leaf == tree_semi_leaves.end()) // if it was not a semi-leaf, it must become one 
                tree_semi_leaves.push_back(parent_node);
        }

        // no need for bound update; but we need to remove the node from leaf nodes and altogether
        tree_leaves.erase(find(tree_leaves.begin(), tree_leaves.end(), leaf_node));
        
        delete leaf_node->corrBCC;
        delete leaf_node;

        return;
    }

    running_bound-=leaf_to_root_bound(leaf_node); // remove contribution of leaf node right away
    if(DEBUG) cout << "Running bound decreased by " << leaf_to_root_bound(leaf_node) << " (contribution of leaf node)" << endl;

    // check if the node is indeed a leaf
    if(!B->isLeaf || leaf_node->children.size() > 0){
        // printMPtree(mptree);
        // printLeaves();
        // printSemiLeaves();
        // cout << "Trying to expand component with ID " << B->myID << ", which has isLeaf="<< B->isLeaf << " and number of children " << leaf_node->children.size() << endl;
        throw invalid_argument("Trying to expand a non-leaf node of the multi-source paths tree (BCC is non leaf or tree node has children)");
    }

    // choose a source at random 
    int sIndex = rand() % B->sources.size();
    pair<int,int> sBpair = B->sources[sIndex];

    // remember multiplicity of s, as it will be the starting one of all its neighbors
    int og_multiplicity = sBpair.second;
    int sB = sBpair.first;

    if(og_multiplicity == -1){ // the BCC had children BCCs
        // printMPtree(mptree);
        // printLeaves();
        // printSemiLeaves();
        // cout << "Trying to expand component with ID " << B->myID << endl;
        throw logic_error("Trying to expand a non-leaf node of the multi-source paths tree (source is art point)");
    }


    // we need to consider two cases according to the size of B->sources!
    // 1) if it is 1, we remove B and substitute it with the new tree of components after the removal of s
    // 2) if it is >1, B stays and we modify a copy newB, which will be a new child of the parent of leaf_node
    // the same for the tree node
    // in any case, the selected source must be removed from the graph's sources
    // if(sBpair.second > 1){
    //     // B->sources[sIndex].second--; // NO:  WE REMOVE THE SOURCE ALTOGETHER!
    //     B->sources.erase(B->sources.begin() + sIndex);
    //     B->bound_multiplier-=og_multiplicity;
    // }
    // else
    //     B->sources.erase(B->sources.begin() + sIndex);

    bool is_copy = false;
    BCC* newB = B; // start with a pointer to B
    mptnode* node_to_replace = leaf_node;
    if(B->sources.size()>1){ // if the number of sources is greater than one, create an actual copy
        newB = new BCC;
        // newB->sources = B->sources; // NO: IT MUST ONLY HAVE THE SOURCE I AM REMOVING AS SOURCE
        newB->sources = {sBpair};
        newB->target = B->target;
        newB->personalBound = B->personalBound;
        newB->prodAncestors = B->prodAncestors;
        newB->bound_multiplier = sBpair.second; // only source is sB
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
        
        is_copy = true;
        B->sources.erase(B->sources.begin() + sIndex); // as for old leaf, remove s as a source and remove its contribution to the multiplier
        B->bound_multiplier-=og_multiplicity;

        // also need to add the new tree node it as child of leaf_node parent
        node_to_replace->parent->children.push_back(node_to_replace);
    }
    // now newB, node_to_replace hold the correct graph to modify

    // if we do not create a copy, remove from leaves
    if(!is_copy)
        tree_leaves.erase(find(tree_leaves.begin(), tree_leaves.end(), leaf_node));

    // if(is_copy)
    //     cout << "Creating a copy" << endl <<flush;

    // remove s from nodes, incident edges of s from edges
    newB->nodes.erase(find(newB->nodes.begin(), newB->nodes.end(), sB));
    vector<int> source_neighbors = newB->edges[sB]; // save neighbors of s for later

    // if(DEBUG){
        // cout << "Neighbors of the source are: "; 
        // for(auto sneigh:source_neighbors)
        //     cout << sneigh << " ";
        // cout << endl << flush;
    // }

    newB->edges.erase(newB->edges.find(sB));
    for (auto neigh: source_neighbors) // also delete from neighbors' list
        newB->edges[neigh].erase(find(newB->edges[neigh].begin(), newB->edges[neigh].end(), sB));

    if(DEBUG) {
        cout << "Erased source "<< sB << "; current graph is " << endl;
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

    // NOTE: IF THE TARGET IS A NEIGHBOR OF THE SOURCE, WE NEED TO INCREASE ITS MULTIPLICITY AS SOURCE OF THE PARENT BCC
    auto neigh_it = find(source_neighbors.begin(), source_neighbors.end(), newB->target);
    if(neigh_it != source_neighbors.end()){ 
        if(DEBUG) cout << "Source of current leaf has target as a neighbor" << endl;
        BCC* source_to_update = leaf_node->parent->corrBCC; // go to BCC of the parent node in the mptree

        bool changed_semi_leaf = false;
        int semi_leaf_index = -1;   
        auto parent_semi_leaf = find(tree_semi_leaves.begin(), tree_semi_leaves.end(), leaf_node->parent); // NOTE: CAN HAPPEN FOR AT MOST ONE SEMI LEAF
        if(parent_semi_leaf != tree_semi_leaves.end()){
            changed_semi_leaf = true; 
            semi_leaf_index = parent_semi_leaf - tree_semi_leaves.begin(); // we indeed changed the contribution of a semi leaf, specifically the one at position semi_leaf_index 
            running_bound -= leaf_to_root_bound(tree_semi_leaves[semi_leaf_index]); // remove from bound before changing multiplicity
            
            if(DEBUG) cout << "Running bound decreased by " << leaf_to_root_bound(tree_semi_leaves[semi_leaf_index]) << " (contribution of parent semi leaf node)" << endl;
        }

        auto update_it = find_if(source_to_update->sources.begin(), source_to_update->sources.end(), [& newB](const pair<int,int> p){return p.first == newB->target;}); // find position of source to update
        // now, check if pair has -1 as second element, change it to og multiplicity, otherwise increase it by og multiplicity. Also, increase multiplier by og multiplicity in any case
        if(source_to_update->sources[update_it - source_to_update->sources.begin()].second == -1)
            source_to_update->sources[update_it - source_to_update->sources.begin()].second = og_multiplicity;
        else
            source_to_update->sources[update_it - source_to_update->sources.begin()].second+=og_multiplicity;

        source_to_update->bound_multiplier+=og_multiplicity;

        if(changed_semi_leaf){
            running_bound+=leaf_to_root_bound(tree_semi_leaves[semi_leaf_index]); // re-add contribution after changing multiplicity
            
            if(DEBUG) cout << "Running bound increased by " << leaf_to_root_bound(tree_semi_leaves[semi_leaf_index]) << " (contribution of parent semi leaf node after adjusting parameters)" << endl;
        }
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
    uint64_t before_bcc = timeMs();
    findBCCs(newB->target, newB, decomposed_BCC, source_neighbors, og_multiplicity);
    

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
        unordered_set<int> curr_nodes = {};

        // for the current BCC we need to unstack until the stack is empty
        // since stack is actually a vector, go through the whole vector
        for(auto stack_el : tree_stack){
            current_BCC->nodes.push_back(stack_el);
            curr_nodes.insert(stack_el);
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
                // if(find(current_BCC->nodes.begin(), current_BCC->nodes.end(), y) != current_BCC->nodes.end()){
                //     current_BCC->edges[x].push_back(y); // add y to edges' vector of x
                //     current_BCC->num_edges++;
                // }

                if(curr_nodes.find(y)!= curr_nodes.end()){
                    current_BCC->edges[x].push_back(y); // add y to edges' vector of x
                    current_BCC->num_edges++;
                }
            }
        }

        current_BCC->num_edges = current_BCC->num_edges/2;

        // PERSONAL BOUND MUST BE CYCLOMATIC MULTIPLIED FOR THE NUMBER OF "VALID" SOURCES = NOT CONTAINING -1
        // int cyclomatic_bound = current_BCC->num_edges - current_BCC->nodes.size() + 1; // initialize personal bound
        current_BCC->personalBound = current_BCC->num_edges - current_BCC->nodes.size() + 1; 
        if(current_BCC->nodes.size() == 2) current_BCC->personalBound = 1;

        bool only_art_pts = true;
        
        if(DEBUG) cout << "Number of nodes in this BCC is " << current_BCC->nodes.size()<<endl;
        // current_BCC->personalBound = cyclomatic_bound;
        for (auto spair : current_BCC->sources)
        {
            if(spair.second != -1){
                only_art_pts = false;
                // current_BCC->personalBound += cyclomatic_bound*spair.second;
                current_BCC->bound_multiplier += spair.second;
            }
        }

        // if(only_art_pts) current_BCC->personalBound = cyclomatic_bound;
        

        // also, set up if it is a leaf: that is, if it ONLY has sources with multiplicity > 0, that is, they are neighbors of the source
        if(find_if(current_BCC->sources.begin(), current_BCC->sources.end(), [&current_BCC](const pair<int, int>& p){return p.second == -1;}) == current_BCC->sources.end()){
            current_BCC->isLeaf = true;
        }
        
        // add the BCC to the final vector
        decomposed_BCC.push_back(current_BCC);
    }

    
    if(DEBUG){
        cout << "BCC found are: "<<  endl;
        for (int i = 0; i < decomposed_BCC.size(); i++){
            cout << "B_"  << i << ": ";
            printBCC(decomposed_BCC[i]);
        }
    }

    // cout << "TIME FOR BCC: " << timeMs() - before_bcc <<endl<<flush;

    // need to build tree and update prodAncestors fields 
    mptnode* new_subtree; // this is a pointer to the root of the new tree; the root corresponds to the only component holding newB->target as target
    int og_leaves = tree_leaves.size(); // save for later vector<mptnode*> root_nodes; // roots of the new subtree
    int og_semi_leaves = tree_semi_leaves.size();
    before_bcc = timeMs();
    update_tree(newB, decomposed_BCC, node_to_replace);
    // cout << "TIME FOR TREE UPDATE: " << timeMs() - before_bcc <<endl<<flush;


    // UPDATE RUNNING BOUND: PROBLEM AS WHEN WE HAVE BOTH AN ART POINT AND A SOURCE IT DOES NOT COUNT AS LEAF   
    //============================= REMOVAL OF LEAF  TO ROOT PATH OF LEAF_NODE PERFORMED AT THE BEGINNING ========================
    // 1) remove leaf-to-root contribution of node_to_replace
    // running_bound -= leaf_to_root_bound(node_to_replace); // note: at beginning, running is 1 and so is the leaf to root of the single node

    // // first, remove B from tree_leaves and from running bound ONLY IF IT IS NOT A COPY
    // if(!is_copy){
    //     tree_leaves.erase(find(tree_leaves.begin(), tree_leaves.end(), leaf_node));
    //     running_bound -= leaf_to_root_bound(leaf_node);
    // }

    // if(is_copy){ // if it is a copy, we also need to recompute the bound for the other copy, as the multiplier has changed
    //     // running_bound -= leaf_to_root_bound(leaf_node); // thus, we subtract one more and add again the contribution of the copy ALREADY DONE BEFORE
    //     // B->bound_multiplier-= og_multiplicity; // also need to decrease the bound multiplier by the multiplicity of the source we removed DONE BEFORE
    //     running_bound+= leaf_to_root_bound(leaf_node); // re-add after multiplicity of source has been decreased 
    // }


    // if(is_copy){ // if it is a copy, we also need to recompute the bound for the other copy, as the multiplier has changed
    //     // running_bound -= leaf_to_root_bound(leaf_node); // thus, we subtract one more and add again the contribution of the copy ALREADY DONE BEFORE
    //     // B->bound_multiplier-= og_multiplicity; // also need to decrease the bound multiplier by the multiplicity of the source we removed DONE BEFORE
    //     running_bound+= leaf_to_root_bound(leaf_node); // re-add after multiplicity of source has been decreased 
    // }

    // if it was a copy, re-add contribution of original node (it will have changed because of source removal)
    if(is_copy){
        running_bound+=leaf_to_root_bound(leaf_node);
        
        if(DEBUG) cout << "Running bound increased by " << leaf_to_root_bound(leaf_node) << " (contribution of leaf node)" << endl; 
    }

    // 2) add contributions of all leaves to it 
    for (auto i = og_leaves; i < tree_leaves.size(); i++){ // only iterates on the new leaves = the ones of the new subtree
        running_bound += leaf_to_root_bound(tree_leaves[i]);
        if(DEBUG) cout << "Running bound increased by " << leaf_to_root_bound(tree_leaves[i]) << " (contribution of new leaf node)" << endl;
    }
    
    for (auto i = og_semi_leaves; i < tree_semi_leaves.size(); i++){ // only iterates on the new semi leaves = the ones of the new subtree
        running_bound += leaf_to_root_bound(tree_semi_leaves[i]);
        if(DEBUG) cout << "Running bound increased by " << leaf_to_root_bound(tree_semi_leaves[i]) << " (contribution of new semi leaf node)" << endl;
    }

    // NOTE: we could have increased the multiplicity of a semi leaf if the source was neighboring the target. In this case, subtract old and re-add new contribution
    // ----- done right away 
    
    // remove the full node and its BCC
    delete node_to_replace->corrBCC;
    delete node_to_replace;

    // if(DEBUG) cout << "Exited tree update procedure" <<endl<<flush;
    // if(DEBUG){
    //     cout << "After updating the tree the new BCCs are: "<<  endl;
    //     for (int i = 0; i < decomposed_BCC.size(); i++){
    //         cout << "B_"  << i << ": ";
    //         printBCC(decomposed_BCC[i]);
    //     }
    // }

    return;
}

long threshold= 1000;

bool assess_paths(mptnode* current_node){
    auto leaf_it = find(tree_leaves.begin(), tree_leaves.end(), current_node);
    int leaf_index = leaf_it - tree_leaves.begin();

    uint64_t assess_start = timeMs();

    while (running_bound < z){

        if(tree_leaves.size() == 0)
            throw logic_error("Leaves' array cannot be empty before terminating");

        if(MAX_TIME>0 && timeMs() - assess_start > MAX_TIME)
            return running_bound >= z;

        // if we are at the root of the tree, we are done
        if(current_node == mptree){
            if(tree_leaves.size() == 1 && current_node->children.size() == 0){
                running_bound = current_node->corrBCC->bound_multiplier;
                return running_bound >= z;
            }
            else{ // pick another leaf
                // this should be impossible; if the root is a leaf, then no more nodes must exist
                throw logic_error("Root added as leaf when tree is non-empty");
                while (current_node == mptree)
                {
                    // choose the next current node at random from the leaves
                    current_node = pick_random_leaf();
                    // leaf_index = rand() % tree_leaves.size();
                    // current_node = tree_leaves[leaf_index];
                }
            }
        }

        if(DEBUG) {
        cout << "Considering component with ID="<< current_node->corrID << endl<< flush; 
        cout << "Component has size " << current_node->corrBCC->nodes.size() <<endl<<flush;
        }


        if(DEBUG){
            printLeaves();
            printSemiLeaves();
        }

        calls_performed++;

        explode(current_node); // explode also removes from leaves
        // cout << endl << endl<< "Finished one explosion"<<endl<<flush;
        // printMPtree(mptree);
        // cout <<flush;

        if(DEBUG){
            printMPtree(mptree);
            cout << endl;
        }

        // choose the next current node at random from the leaves
        current_node = pick_random_leaf();
        // leaf_index = rand() % tree_leaves.size();
        // current_node = tree_leaves[leaf_index];
        // tree_leaves.erase(tree_leaves.begin() + leaf_index);
 

        // cout << "Updated running bound: "<< running_bound <<endl;

        
        if(running_bound>threshold){
            cout << "lb="<< running_bound<< "; calls=" << calls_performed << "; t="<< timeMs() - assess_start << "\t\t";
            threshold*=10;
        }
    }

    return true; // exited while bound < z
}

int main(int argc, char** argv) { 
    if(argc < 5){
        cout << "USAGE: " << argv[0] << " <graph-filename> <source> <target> <z> [MAX_TIME] [OUTPUT_FILENAME]\n";
        return 0;
    }

    string outname;


    if(argc >= 6) MAX_TIME = atoi(argv[5])*1000;
    

    s = atoi(argv[2]);
    t = atoi(argv[3]);
    z = atoi(argv[4]);

    int t_og = t;

    if(s==t)
        throw invalid_argument("Source and target must be different");
    if(z<0)
        throw invalid_argument("Parameter z must be a positive integer");

    create_graph(argv[1]);
    
    #ifdef RANDOM_NODE_EXTRACTION
        srand(time(NULL));
    #endif

    // srand(13);

    mptnode* graph_node = mptree->children[0];

    if(DEBUG){
        cout << "Original graph: " << endl<<flush;
        for(auto x : graph_node->corrBCC->nodes){
            cout << x << ": ";
            for(auto y : graph_node->corrBCC->edges[x]){
                cout << y << " ";
            }
            cout << endl;
        }
    }

    // cout << "Starting running bound: "<< running_bound << endl;
    // explode(graph_node);
    
    // printMPtree(mptree);
    // cout << endl;

    // cout << "Ending running bound: "<< running_bound <<endl;
    tree_leaves.push_back(graph_node);

    uint64_t start_time = timeMs();
    bool enough = assess_paths(graph_node);
    uint64_t duration = (timeMs() - start_time);

    // cout << "Elapsed time: " << duration << "ms"  << endl;

    // if(enough)
    //     cout << "The paths are at least z=" << z << endl;
    // else    
    //     cout << "The paths are exactly " << running_bound << " < z=" << z << endl;

    if(argc >= 7){
        ofstream output_graph;
        output_graph.open(argv[6]);
        output_graph << endl << "Assessment: "<< argv[1] << " "<< N << " " << M <<  " " << s << " " << t_og << "; " << duration << " " << calls_performed  << " " << running_bound << " "<< z << " " << running_bound/duration << endl;
        output_graph.close();
    }
    else
        cout << endl << "Assessment: "<< argv[1] << " "<< N << " " << M <<  " " << s << " " << t_og << "; " << duration << " " << calls_performed  << " " << running_bound << " "<< z << " " << running_bound/duration <<endl;


    delete_all(mptree);

    return 0;
}