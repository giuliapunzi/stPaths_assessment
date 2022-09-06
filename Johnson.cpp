#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstdint>

using namespace std;

////////////////// JOHNSON STUFF //////////////////////////

#include <stack>
#include <unordered_set>

// #define DEBUG true
#define DEBUG false

vector<int> blocked;
vector<unordered_set<int>*> B;
long all_paths = 0;
long MAX_TIME = -1;

int SOURCE;
int TARGET;

///////////////////////////////////////////////////////////

vector<int> deleted;
vector<int> current_degree;
typedef vector<vector<int>> graph;
graph G;


bool is_edge(int u, int v);


uint64_t timeMs() {
  using namespace std::chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}


// create graph from file filename (WE REMOVE MULTI-EDGES AND SELF LOOPS)
// NDE format
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
    // reachable.resize(N, 1);

    // initialize reachable vector to be all 1 and degree to be zero
    for(int i = 0; i < G.size() ; i++){
        current_degree[i] = 0;
        // reachable[i]=1;
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
 

// create graph from file filename 
// WE REMOVE SELF LOOPS
graph create_graph_aa(char* filename)
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
        }
        
    }

    cout << "Input graph has " << N << " nodes and " << real_edges << " edges. "<< endl;

    fclose(input_graph);
    return G;
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


// returns degree of node u
// remember: need to take into account deleted edges
inline int degree(int u)
{
    int deg = 0;
    if(deleted[u])
        return -1;
    else{
        for(int i = 0; i< G[u].size(); i++){
            if(!deleted[G[u][i]])
                deg++;
        }
    }

    return deg;
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

inline void printDeleted(){
    for(int i= 0; i<deleted.size(); i++){
        if(deleted[i])
            cout << i << " ";
    }
    cout << endl;
}

// ++++++++++++++++++++++++++++++++++ GRAPH MODIFIERS ++++++++++++++++++++++++++++++++

// remove node u from the graph, by marking deleted vector to 0
inline void remove_node(int u)
{
    deleted[u] = 1;
    return;
}


// ++++++++++++++++++++++++++++++++ GRAPH VISITS +++++++++++++++++++++++++++++++

// recursive DFS procedure from node u
void DFS(int u, vector<bool> &visited){
    // cout << "Entering DFS for " << u << endl << flush;
    visited[u] = true;

    for(int i = 0; i < G[u].size(); i++)
    {
        if(!visited[G[u][i]] && !deleted[G[u][i]])
            DFS(G[u][i], visited);
    }

    return;
}

int visits_performed;

// starts a visit from t, and marks as good neighbors the neighbors of s that
// are reached through the visit. Outputs the vector of these good neighbors.
vector<bool> check_neighbors(int s, int t){
    vector<bool> visited(G.size());
    visits_performed++;

    // DELETION OF S PERFORMED BEFORE CALL TO FUNCTION
    // before starting, mark s as deleted
    // deleted[s] = 1;

    // initialize all deleted nodes as visited
    for(int i = 0; i< visited.size(); i++){
        if(deleted[i])
            visited[i] = true;
        else
            visited[i] = false;
    }

    // for(auto x : visited)
    //     cout << x << " " << flush;
    // cout << endl;

    // launch DFS from node t
    DFS(t, visited);

    vector<int> neigh = neighbors(s);
    // find out which neighbors of s have been visited, and output them
    vector<bool> good_neighbors(neigh.size());
    for(int i = 0; i < neigh.size(); i++){
        if(visited[neigh[i]])
            good_neighbors[i] = true;
            // good_neighbors.push_back(G[s][i]);
        else
            good_neighbors[i] = false;
            
    }
    
    // cout << "Before outputting good neighbors: size is " << good_neighbors.size() << ", while size of neighbors is "<< neigh.size() << endl<< flush;

    // undo deletion of s 
    // deleted[s] = 0;

    return good_neighbors;
}

// global variable used to count the number of paths
// also count the total length of the paths up to now 
// (can be substituted with full enumeration)
unsigned long count_paths;
unsigned long long total_length;
int curr_path_len;
unsigned long good_diff_len;

// global variable used to find how many dead ends there are
unsigned long dead_ends;
unsigned long long dead_total_len; // total length of dead ends
unsigned long dead_diff_len; // edges only belonging to dead ends; increase by 1 every time we backtrack

// constant which limits the number of function calls
// plus global variable that takes into account the number of calls performed
const long MAX_CALLS = 500000000000;
// const long MAX_TIME = 60000; 
int calls_performed;
uint64_t start_time;



void enumerate_paths(int s, int t){
    count_paths = 0;
    total_length = 0;
    dead_ends = 0;
    calls_performed = 0;
    curr_path_len = -1;
    good_diff_len = 0;
    dead_diff_len = 0;
    dead_total_len = 0;
    visits_performed=0;
    // paths(s,t);
    good_diff_len--; // source returned true and thus added one 

    return;
}





// void badUNBLOCK (int u){
//         blocked[u] = 0;
//         while(! B[u].empty()){
//                 int w = B[u].pop(); // delete w from B(u);
//                 if(blocked[w]) UNBLOCK(w);
//             }
//     }


// assumes B[u] is a set, not recursive and removal safe
void UNBLOCK(int u){
        stack<int> todo;
        todo.push(u);
        int x;

        if(DEBUG) cout << "UNBLOCK----------\n";
        
        while(!todo.empty()){
            x = todo.top();
            todo.pop();
            if(DEBUG) cout << x << " (";
            blocked[x] = 0;
            for(int y : *(B[x])){
                if(blocked[y]){
                    todo.push(y);
                    blocked[y] = 0;
                    if(DEBUG) cout << y << " ";
                }
            }
            if(DEBUG) cout << ")\n";
            B[x]->clear();
        }

        if(DEBUG) cout << "END-----UNBLOCK\n";
    }

long time_evals = 0;
long eval_resolution = 10000;
bool abort_alg = false;

short PATHS(int v){

    if(abort_alg) return true;
    else if (MAX_TIME > 0 && time_evals%eval_resolution == 0){
        if(timeMs() - start_time >= MAX_TIME){
            abort_alg=true;
            return true;
        }
    }
    time_evals++;

    // if(MAX_TIME > 0 && (timeMs() - start_time > MAX_TIME)) return 0;

    short f = 0;
//      stack v; //enqueue v to path
        blocked[v]= 1;
//L1:
        for(auto w : G[v]){
            if(w == TARGET) {
//                  output circuit composed of stack followed by TARGET; // commented
                    all_paths++;
                    f = 1;
            }
            else if(!blocked[w]){
                if(PATHS(w)) f=1;
            }
        }
//L2:
        if(f){
            UNBLOCK(v);
        } else{
            for(auto w : G[v]){
                B[w]->insert(v); // maybe check first? //if v not in B(w) then put v on B(w);
            }
        }
//      unstack v; // dequeue v from solution
        return f;
}



int main(int argc, char * argv[]){
    
    if(argc < 4){
        cout << "USAGE: ./" << argv[0] << " <graph> <source> <targer> [MAX_TIME]\n";
        return 0;
    }

    if(argc >= 5) MAX_TIME = atoi(argv[4]);

    SOURCE = atoi(argv[2]);
    TARGET = atoi(argv[3]);

    cout << "G: " << argv[1] << " S: " << SOURCE << " T: " << TARGET;

    if(MAX_TIME > 0) cout << " MAX_TIME: " << MAX_TIME << "ms";
    cout << endl;

    char* input_filename = argv[1];
    create_graph(input_filename);

    // printGraph();


    // find max degree of graph
    int maxdeg = 0;
    for(int u=0; u < G.size(); u++){
        if(maxdeg< degree(u)){
            maxdeg = degree(u);
        }
    }

    // here we also find the number of edges
    int numedges = 0;
    int numnodes = G.size();

    for(int node = 0; node < G.size(); node++){
        numedges += G[node].size();
    }
    numedges = numedges/2;

    cout << "Graph has maximum degree " << maxdeg << " and n=" << G.size() << endl; 

    start_time = timeMs();
    // standard: s = 0, t=last node

    int n = G.size();
    for(int i = 0; i < n; i++){
        blocked.push_back(0);
        B.push_back(new unordered_set<int>());
    }
    PATHS(SOURCE);

    uint64_t duration = (timeMs() - start_time);

    cout << "DONE! time: " << duration << "ms, paths: " << all_paths << endl;

/*****
    enumerate_paths(0, G.size()-1);
    uint64_t duration = (timeMs() - start_time);

    cout << endl <<  "Elapsed time: " << duration << " sec; calls performed are " << calls_performed << endl;
    cout << "Visits performed are  " << visits_performed << endl;
    cout << "Paths found are " <<count_paths << " ; their total length is "<< total_length << " and their partial length is " << good_diff_len << endl;
    cout << "Dead ends are " << dead_ends << " ; their total length is " << dead_total_len << " and their partial length is " << dead_diff_len <<endl;

    // reporting to file
    ofstream output_file; 
    output_file.open("output-v0.txt", ios::app);
    output_file << "-----------------------------------------------------"<< endl;
    output_file << "Output for graph with " << numnodes << " nodes, " << numedges << " edges and max degree " << maxdeg << " (" << input_filename << ")"<< endl;
    output_file << calls_performed << " calls performed in " << duration << " secs (MAX_CALLS = " << MAX_CALLS << ")" << endl;
    output_file << "Visits of the graph performed are  " << visits_performed << endl;
    output_file << "Paths found are " <<count_paths << " for a total length of " << total_length << " and a partial length of " << good_diff_len << endl;
    output_file<< "Dead ends are " << dead_ends << " for a total length of "<< dead_total_len << " and a partial length of " << dead_diff_len << endl;
    output_file << "-----------------------------------------------------"<< endl<<endl<<endl;
    output_file.close();


    // reporting to file for exhaustive
    // output_file.open("exhaustive-comparison.txt", ios::app);
    // output_file << "------------------------VERSION 0-----------------------------"<< endl;
    // output_file << "Output for graph with " << numnodes << " nodes, " << numedges << " edges and max degree " << maxdeg << " (" << input_filename << ")"<< endl;
    // output_file << calls_performed << " calls performed in " << duration << " secs (MAX_CALLS = " << MAX_CALLS << ")" << endl;
    // output_file << "Visits of the graph performed are  " << visits_performed << endl;
    // output_file << "Paths found are " <<count_paths << " for a total length of " << total_length << " and a partial length of " << good_diff_len << endl;
    // output_file<< "Dead ends are " << dead_ends << " for a total length of "<< dead_total_len << " and a partial length of " << dead_diff_len << endl;
    // output_file << "-----------------------------------------------------"<< endl<<endl<<endl;
    // output_file.close();
******/

    return 0;
}



// integer list array A(n), B(n); logical array blocked (n); integer s;
//     empty stack;
//     s:=l;
//     while s < n do
//         begin
//             A:= adjacency structure of strong component K with least
//             vertex in subgraph of G induced by {s, s+ 1, n};
//             if A then
//                 begin
//                     s := least vertex in V;
//                     for iVu, do
//                         begin
//                             blocked(i) :-- false;
//                             B(i) :-- ,;
//                         end;
// //L3:
//                     dummy :- CIRCUIT(s);
//                     s:=s+l;
//                 end
//             else s n;
//         end
