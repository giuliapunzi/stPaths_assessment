#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stack>
#include <stdexcept>
#include <chrono>
#include <cstdint>
#include <stack>
#include <sys/resource.h>

using namespace std;

vector<int> deleted;
vector<int> current_degree; // keep global degree vector, updating it as deletions and insertions go

stack<int> del_stack; // stack of nodes that have been deleted
vector<vector<int>> G;

vector<bool> reachable; // marks nodes that have been visited by the DFS
long visits_performed;

// global variable used to count the number of paths
// also count the total length of the paths up to now 
// (can be substituted with full enumeration)
long long count_paths;
long long total_length;
int curr_path_len;
long long good_diff_len;
long deleted_w_caterpillar;
long time_reachability;
long time_caterpillar;

long MAX_TIME=0;

uint64_t start_time;
long long calls_performed;

uint64_t z;


// vector<vector<int>> BCC_stack;
// vector<int> target_stack;

bool is_edge(int u, int v);



long time_evals = 0;
long eval_resolution = 1000;
bool abort_alg = false;


uint64_t timeMs() {
  using namespace std::chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}


// create graph from file filename (WE REMOVE MULTI-EDGES AND SELF LOOPS)
void create_graph(char* filename)
{
    FILE* input_graph = fopen(filename, "r");
    int N;

    // file contains number of nodes, number of edges at fist line
    // and then one edge per line
    fscanf(input_graph, "%d", &N);

    int real_edges = 0;

    G.resize(N);
    deleted.resize(N);
    current_degree.resize(N, 0);
    reachable.resize(N, 1);

    // initialize reachable vector to be all 1 and degree to be zero
    for(int i = 0; i < G.size() ; i++){
        reachable[i] = 1;
        current_degree[i] = 0;
    }
    
    int u, v;
    // we need to skip the first N rows
    for(int i = 0; i < N; i++)
        fscanf(input_graph, "%d %d", &u, &v);
    
    // for(int i=0; i<M; i++)
    while(fscanf(input_graph,"%d %d",&u, &v) != EOF)
    {
        // make sure no self-loops or multiedges are created 
        if (u != v && !is_edge(u,v)){
            G[u].push_back(v);
            G[v].push_back(u);
            real_edges++;
            current_degree[u]++;
            current_degree[v]++;
        }
        
    }

    // cout << "Input graph has " << N << " nodes and " << real_edges << " edges. "<< endl;


    fclose(input_graph);
    return;
}
 

// create graph from file filename FOR OLDER GRAPH FORMAT
void create_graph_old(char* filename)
{
    FILE* input_graph = fopen(filename, "r");
    // input_graph.open(filename);
    
    int N, M;

    // file contains number of nodes, number of edges at fist line
    // and then one edge per line
    // input_graph >> N >> M;
    fscanf(input_graph, "%d %d", &N, &M);

    int real_edges = 0;

    G.resize(N);
    deleted.resize(N);
    current_degree.resize(N, 0);
    reachable.resize(N, 1);

    // initialize reachable vector to be all 1 and degree to be zero
    for(int i = 0; i < G.size() ; i++){
        reachable[i] = 1;
        current_degree[i] = 0;
    }
        

    int u, v;
    for(int i=0; i<M; i++)
    {
        fscanf(input_graph, "%d,%d", &u, &v);
        // input_graph >> u >> v;
        // make sure no self-loops or multiedges are created 
        if (u != v && !is_edge(u,v)){
            G[u].push_back(v);
            G[v].push_back(u);
            real_edges++;
            current_degree[u]++;
            current_degree[v]++;
        }
        
    }

    cout << "Input graph has " << N << " nodes and " << real_edges << " edges. "<< endl;


    fclose(input_graph);
    return;
}


// ++++++++++++++++++++++++++++++++++ GRAPH REPORTING ++++++++++++++++++++++++++++++++

// checks if a given edge (u,v) belongs to graph G
inline bool is_edge(int u, int v)
{   
    if(deleted[u] || deleted[v])
        return false;
    
    if(find(G[u].begin(), G[u].end(), v) == G[u].end())
        return false;
    else
        return true;
}


inline bool is_neighbor(int u, int v){
    if(deleted[v])
        return false;
    
    return find(G[u].begin(), G[u].end(), v) != G[u].end();
}


inline int degree(int u){
    return current_degree[u];
}

// outputs the vector of (non removed) neighbors of u
// NOTE: u can be a removed node
inline vector<int> neighbors(int u)
{
    vector<int> neigh; 
    for(int i = 0; i<G[u].size(); i++){
        if(!deleted[G[u][i]])
            neigh.push_back(G[u][i]);
    }

    return neigh;
}


// print graph as list of adjacency lists
inline void printGraph()
{
    cout << endl;
    // for each node, look through its adjacency list
    for(int i =0; i < G.size(); i++)
    {
        if(!deleted[i]){
            cout << i << ": ";
            for (int j = 0; j< G[i].size(); j++){
                if(!deleted[G[i][j]])
                    cout << '(' << i << ", " << G[i][j] << "); ";
            }
            cout << endl;
        }
        
    }
    cout << endl;
}


// ++++++++++++++++++++++++++++++++++ GRAPH MODIFIERS ++++++++++++++++++++++++++++++++

// remove node u from the graph, by marking deleted vector to 0
inline void remove_node(int u)
{
    deleted[u] = 1;
    del_stack.push(u);

    for (auto v : G[u])
        current_degree[v]--;
    
    return;
}

// reinserts top of the stack node in the graph
inline void reinsert_node()
{
    int u = del_stack.top();
    del_stack.pop();
    deleted[u] = 0;

    for (auto v : G[u])
        current_degree[v]++;
    
    return;
}



inline void remove_simple(int u){
    deleted[u] = 1;

    for (auto v : G[u])
        current_degree[v]--;
    
    // cout << "Removing simple " << u << endl;
    return;
}

inline void reinsert_simple(int u){
    deleted[u] = 0;

    for (auto v : G[u])
        current_degree[v]++;
    
    // cout << "Reinserting simple " << u << endl;
    return;
}


// ++++++++++++++++++++++++++++++++ GRAPH VISITS +++++++++++++++++++++++++++++++


vector<bool> visited;
vector<int> disc;
vector<int> low;
vector<int> parent; 
vector<int> cat_stack;
int visit_time;
int current_s;
int last_art;
bool found_s;
bool root_art_pt;

void find_artpts(int s, int u)
{
    // cout << "Inserting "<< u << " in cat stack"<< endl<<flush;
    cat_stack.push_back(u);
    // Count of children in DFS Tree
    int children = 0;
 
    // Mark the current node as visited
    visited[u] = true;

    if(u== s)
        found_s = true;
 
    // Initialize discovery time and lowpoint value
    disc[u] = low[u] = ++visit_time;
    
    bool root_found = false; // becomes true when the root finds s
    bool good_for_current_BCC;    

    // Go through all non-deleted neighbors of u
    for (auto v : G[u]) {
        // cout <<v << " is a neighbor of " << u << endl;
        if(!deleted[v]){
            // If v is not visited yet, then make it a child of u
            // in DFS tree and recur for it
            if (!visited[v]) {
                good_for_current_BCC = false;
                parent[v] = u;
                children++;
                find_artpts(s, v);

                // if we are the root and we just found s, v is the only good neighbor
                if(parent[u] == -1 && found_s){
                    // if(!root_found && found_s){ // the nodes in the stack right now are the 'correct BCC' for t
                    // BCC_stack.push_back({v});

                    // it could be the second BCC for t so we need to unstack until v
                    int x = cat_stack.back();
                    while(x != v){ // If we are in the "right" BCC add it to the vector   
                        // BCC_stack.back().push_back(x);
                        cat_stack.pop_back();
                        x = cat_stack.back();
                    }
                    cat_stack.pop_back(); 
                    // int ind = cat_stack.size() -1;
                    // while(cat_stack[ind] != v)
                    // {   
                    //     // if the element of the stack is a neighbor of t 
                    //     // if (find(G[u].begin(), G[u].end(), cat_stack[ind])!= G[u].end())
                    //     if (is_neighbor(u, cat_stack[ind])){
                    //         root_correct_neigh.push_back(cat_stack[ind]);
                    //     }  
                    //     ind--;
                    // }
                    // root_found = true;
                    
                    // WE SHOULD BE ABLE TO RETURN HERE!
                    return;
                
                }
    
                // Check if the subtree rooted with v has
                // a connection to one of the ancestors of u
                low[u] = min(low[u], low[v]);
    
                // If u is not root and low value of one of
                // its child is more than discovery value of u.
                if (parent[u] != -1 && low[v] >= disc[u]){ // here is where I close my articulation point
                    // cout << "Closed an articulation point " << endl << flush;
                    int check = cat_stack.size()-1;
                    // first, find out if good for current BCC by scanning only the current BCC
                    while(cat_stack[check] != v)
                    {
                        if (cat_stack[check]== current_s){
                            good_for_current_BCC = true;
                        }  
                        check--;
                    }
                    if(v == current_s) 
                        good_for_current_BCC = true;
                    
                    // NEW: SET THE CURRENT ART POINT AS S!
                    if(good_for_current_BCC){
                        // if at this point the current source was s, this is the last art point
                        // NEED TO CHECK THAT IT IS DIFFERENT FROM S
                        // if(current_s == s && u != s)
                        //     last_art = u;
                        
                        // BCC_stack.push_back({v}); // initialize the BCC as v

                        // remove stuff from stack
                        int x = cat_stack.back();
                        while(x != v){ // If we are in the "right" BCC add it to the vector                            
                            // BCC_stack.back().push_back(x);
                            cat_stack.pop_back();
                            x = cat_stack.back();
                        }
                        cat_stack.pop_back(); 

                        // if(!is_neighbor(u, v))
                        //     throw logic_error("v should be a neighbor of u!!!");

                        current_s = u; 

                        // push the articulation point to the stack, and its vector of good neighbors in the corresponding
                        // position of the target_neighbors vector. 
                        // temp_target_stack.push_back(u);
                        // temp_target_neigh.push_back(good_neigh_current_target);
                        // target_stack.push_back(u);
                    }
                        
                    if(!good_for_current_BCC){ // if I am not good, I need to delete my neighbors that have discovery time greater than v
                        // remove stuff and delete
                        int x = cat_stack.back();
                        while(x != v){ // NEW: ONLY DELETE NEIGHBORS, SO MORE CHECKS BUT LESS WORK OVERALL
                            cat_stack.pop_back();
                            if(is_neighbor(u,x)){
                                remove_node(x);
                                // deleted_w_caterpillar++;
                                // cout << "[caterpillar removed " << x << "]" << endl;
                            }
                            
                            x = cat_stack.back();
                        }
                        cat_stack.pop_back();

                        if(!is_neighbor(u, v))
                            throw logic_error("v should be a neighbor of u!!!");

                        remove_node(x); // v is always neighbor
                        // deleted_w_caterpillar++;
                        // cout << "[caterpillar removed " << x << "]" << endl;
                    }
                }
            }
    
            // Update low value of u for parent function calls.
            else if (v != parent[u])
                low[u] = min(low[u], disc[v]);
        }
    }
    
    return;
    // If u is root of DFS tree and has two or more children.
    // In this case, we don't actually care about removing other BCCs, as they will never be explored
    // if (parent[u] == -1 && children > 1){
    //     root_art_pt = true;

    //     // remove neighbors of t that are bad, i.e. the neighbors not in root_correct_neigh
    //     // these can be removed, as they do not change as we go deeper in the recursion, only when we return normally. 
    //     // for(auto neigh : neighbors(u)){
    //     //     // cout << "Considering " << neigh << endl;
    //     //     if(find(root_correct_neigh.begin(), root_correct_neigh.end(), neigh) == root_correct_neigh.end()){ // need to check if neighbor's parent is t
    //     //         // cout << "Removing " << neigh << endl;
    //     //         remove_node(neigh);
    //     //         deleted_w_caterpillar++;
    //     //         // cout << "[caterpillar removed " << neigh << "]" << endl;
    //     //     }
    //     // }
    // }
}

// start visit for finding articulation points from t
void find_caterpillar(int s, int t)
{   
    cat_stack.erase(cat_stack.begin(), cat_stack.end());
    visit_time = 0;
    visits_performed++;
    found_s = false;
    current_s = s;
    last_art = -1;

    for(int i = 0; i < G.size(); i++){
        visited[i] = false;
        parent[i] = -2;
    }
    
    parent[t] = -1;

    // cout << "About to find art points from "<< s << " to "<< t;
    // printGraph();
    uint64_t start = timeMs();

    // sizes of the stacks
    // int size_t_stack = target_stack.size();
    // int size_C_stack = BCC_stack.size();


     // only interested in the ones from s to t = caterpillar
    find_artpts(s, t);

    // after the artpoints call, we added stuff IN REVERSE ORDER from size_C_stack+1 and size_t_stack+1 onwards.    
    // reverse(target_stack.begin() + size_t_stack, target_stack.end());
    // reverse(BCC_stack.begin() + size_C_stack, BCC_stack.end());

    // if(target_stack.size() != BCC_stack.size()-1)
    //     throw logic_error("Wrong number of targets wrt BCCs when exiting caterpillar");


    time_caterpillar += (timeMs() - start);
    return;
}



// recursive DFS procedure from node u
void DFS(int u){
    // cout << "Entering DFS for " << u << endl << flush;
    reachable[u] = 1;

    for(int i = 0; i < G[u].size(); i++)
    {
        // reachable being equal to zero means that the node has not been deleted nor already visited
        if(!reachable[G[u][i]])
            DFS(G[u][i]);
    }

    return;
}




void DFS_iter(int u){
    // initialize reachable
    for(int i = 0; i< reachable.size(); i++)
        reachable[i] = false;

    stack<int> DFS_stack;
    DFS_stack.push(u);

    while (!DFS_stack.empty())
    {
        int v = DFS_stack.top();
        DFS_stack.pop();
        reachable[v] = true;

        for (auto x : G[v])
        {
            if(!reachable[x]) // and not dleted
                DFS_stack.push(x);
        }
        
    }
    
}


// starts a visit from t, and marks as reachable the nodes that
// are reached through the visit.
void reachability_check(int t){
    reachable.resize(G.size());
    visits_performed++;

    // initialize all deleted nodes as visited
    for(int i = 0; i< reachable.size(); i++){
        if(deleted[i])
            reachable[i] = 1;
        else
            reachable[i] = 0;
    }

    uint64_t start = timeMs();
    // launch DFS from node t
    DFS_iter(t);
    time_reachability += (timeMs() - start);

    // go through all nodes of the graph and deleted the ones with reachable value = 0
    // delete means both mark deleted[u] = 0 and add them to the stack of deleted nodes
    for(int u = 0; u < G.size(); u++){
        // here we need the differentiation between -1 and 0: otherwise we add to the stack nodes already removed
        // NO: CHECK IF DELETED BY USING REGULAR VECTOR!
        if(!reachable[u] && !deleted[u]){
            // deleted[u] = true;
            // del_stack.push(u);
            remove_node(u);
        }
    }

    return;
}

uint64_t threshold= 10000;

// paths must return the status, either success or fail
// we do so by returning true/false: true = success
void paths_classic_BCC(int u, int t){
    struct rusage usage_output;
    curr_path_len++;

    if(count_paths >= z)
        return;

    // if(calls_performed >= MAX_CALLS)
    //     return;

    if(MAX_TIME>0 && timeMs() - start_time >= MAX_TIME)
        return;
    
    getrusage(RUSAGE_SELF, &usage_output); // measure memory usage 
    if(usage_output.ru_maxrss>400000000) // if usage greater than about 400GB, return
        return;

    calls_performed++;

    // if(calls_performed % 1000000 == 0)
    //     cout << "*" << flush;

    // base case
    if(u== t){
        count_paths++;
        total_length+=curr_path_len;
        curr_path_len--;
        good_diff_len++;
        return;
    }

    // at every step, we perform the visit which deletes non-reachable nodes
    remove_node(u); // I need to remove the node to ensure correct placement on deleted stack

    reinsert_simple(u);
    find_caterpillar(u, t);
    remove_simple(u);

    if(degree(u) == 0){
        throw logic_error("Error: Found a dead end in correct algorithm");
    }

    // we have non-deleted neighbors; explore them
    // if(degree(u) > 0){ // don't need this anymore
    for(auto v: G[u]){ 
        if(!deleted[v]){ // we take the next non-deleted element of G[u], noting that these deleted elements dynamically change during the for loop
            paths_classic_BCC(v, t);
        }
        
    }

    // rimettere le cose a posto       
    // pop stack until u (included) and mark as not deleted
    while(del_stack.top() != u)
        reinsert_node();
    
    reinsert_node(); // here we are inserting u 


    if(timeMs() - start_time > threshold){
    // if(calls_performed > it_num){ // every 10 iterations
        getrusage(RUSAGE_SELF, &usage_output); // measure memory usage 
        cout << timeMs() - start_time << "\t" << count_paths << "\t"<< calls_performed << "\t" << usage_output.ru_maxrss << endl;

        threshold+=10000;
    }

    curr_path_len--;
    good_diff_len++;

    return;
}

void enumerate_paths_classical(int s, int t){
    count_paths = 0;
    total_length = 0;
    calls_performed = 0;
    curr_path_len = -1;
    good_diff_len = 0;
    visits_performed= 0;
    paths_classic_BCC(s,t);
    good_diff_len--; // source returned true and thus added one 

    return;
}


int main(int argc, char* argv[]){ 
    if(argc < 5){
        cout << "USAGE: " << argv[0] << " <graph-filename> <source> <target> <z> [MAX_TIME] [OUTPUT_FILENAME]\n";
        return 0;
    }

    string outname;


    if(argc >= 6) {
        MAX_TIME = atoi(argv[5])*1000;
    }

    int s = atoi(argv[2]);
    int t = atoi(argv[3]);
    // z = atoi(argv[4]);
    // z= atoll(argv[4]);

    z=stoull(argv[4]);

    if(s==t)
        throw invalid_argument("Source and target must be different");
    if(z<0)
        throw invalid_argument("Parameter z must be a positive integer");

    create_graph(argv[1]); // initialize graph


    disc.resize(G.size());
    low.resize(G.size());
    visited.resize(G.size());
    parent.resize(G.size());

    // find max degree of graph
    // int maxdeg = 0;
    // for(int u=0; u < G.size(); u++){
    //     if(maxdeg< degree(u)){
    //         maxdeg = degree(u);
    //     }
    // }

    // here we also find the number of edges
    int numedges = 0;
    int numnodes = G.size();

    for(int node = 0; node < G.size(); node++){
        numedges += G[node].size();
    }
    numedges = numedges/2;

    // cout << "Graph has maximum degree " << maxdeg << endl; 

    // int s , t;
    // cout << "Insert value for s from 0 to " << numnodes-1 << ": ";
    // cin >> s;
    // cout << "Insert value for t from 0 to " << numnodes-1 << ": ";
    // cin >> t;

    // initialize all nodes as non-reachable
    for(int i = 0; i< reachable.size(); i++)
        reachable[i] = 0;

    // chech reachability of s from t
    DFS_iter(t);

    if(!reachable[s]){
        cout << "Node s is not reachable from t" << endl;
        return 0;
    }

    // mark as deleted the non-reachable nodes
    for(int i =0;i < G.size();i++){
        if(!reachable[i])
            deleted[i]= 1;
    }

    // target_stack.erase(target_stack.begin(), target_stack.end());
    // BCC_stack.erase(BCC_stack.begin(), BCC_stack.end());
    // target_stack.push_back(t);

    deleted_w_caterpillar = 0;
    find_caterpillar(s,t); // last art is now set up
    
    time_reachability = 0;
    time_caterpillar = 0;

    start_time = timeMs();

    enumerate_paths_classical(s, t);
    uint64_t duration = (timeMs() - start_time);

    struct rusage usage_output;
    int usageint = getrusage(RUSAGE_SELF, &usage_output); // measure memory usage 

    if(argc >= 7){
        ofstream output_graph;
        output_graph.open(argv[6]);
        output_graph << endl << "ClassicBCC: "<< argv[1] << " "<< numnodes << " " << numedges <<  " " << s << " " << t << "; " << duration << " " << calls_performed  << " " << count_paths << " " << z << " " << count_paths/duration << endl;
        output_graph << "Memory usage: " << usage_output.ru_maxrss/1000 << "MB"<< endl;
        
        output_graph.close();
    }
    else {
        cout << endl << "ClassicBCC: "<< argv[1] << " "<< numnodes << " " << numedges <<  " " << s << " " << t << "; " << duration << " " << calls_performed  << " " << count_paths << " " << z << " " << count_paths/duration << endl;
        cout << "Memory usage: " << usage_output.ru_maxrss/1000 << "MB"<< endl;
    }


    cout<< endl << endl;

    return 0;
}
