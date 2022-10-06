// extract a subgraph (with corresponding caterpillar and source and target) of at most required size from a file
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <stack>

using namespace std;

// #define DEBUG true
#define DEBUG false

int N,newN;
int MAX_N;
int M = 0;
int newM=0;
int s,t;
string outname;

vector<vector<int>> G;

vector<bool> visited;
vector<int> disc;
vector<int> low;
vector<int> parent; 
vector<int> cat_stack;
vector<int> to_remove;
vector<int> new_names;
vector<int> reachable;

int visit_time;
int current_s;
int last_art;
bool found_s=false;


// create graph from NDE file filename (REMOVE MULTI-EDGES AND SELF LOOPS)
void create_graph(char* filename)
{
    FILE* input_graph = fopen(filename, "r");

    // file contains number of nodes, number of edges at fist line
    // and then one edge per line
    fscanf(input_graph, "%d", &N);

    // if(DEBUG) cout << "Input graph has " << N << " nodes."<< endl;

    if(MAX_N >= N){
        cout << "Input graph has " << N << "<" << MAX_N << "nodes. "<< endl;
        throw invalid_argument("Cannot extract subgraph");
    }

    G.resize(MAX_N);
    visited.resize(MAX_N);
    parent.resize(MAX_N);
    disc.resize(MAX_N);
    low.resize(MAX_N);
    to_remove.resize(MAX_N);
    new_names.resize(MAX_N);
    reachable.resize(MAX_N);
    
    int u, v;
    // we need to skip the first N rows, plus add 0,...,N to nodes
    for(int i = 0; i < N; i++){
        fscanf(input_graph, "%d %d", &u, &v);
    }
    
    // for(int i=0; i<M; i++)
    while(fscanf(input_graph,"%d %d",&u, &v) != EOF)
    {
        // if(DEBUG) cout << "Considering pair " << u << " " << v << endl;
        // make sure no self-loops or multiedges are created 
        if (u != v && u < MAX_N && v < MAX_N && find(G[u].begin(), G[u].end(), v) == G[u].end()){
            G[u].push_back(v);
            G[v].push_back(u);
            M++;
        }   
    }

    // cout << "Input graph has " << N << " nodes and " << M << " edges. "<< endl;

    fclose(input_graph);
    return;
}

// print graph as list of adjacency lists
inline void printGraph()
{
    cout << endl;
    // for each node, look through its adjacency list
    for(int i =0; i < G.size(); i++)
    {
        if(!to_remove[i]){
            cout << i << ": ";
            for (int j = 0; j< G[i].size(); j++){
                if(!to_remove[G[i][j]])
                    cout << '(' << i << ", " << G[i][j] << "); ";
            }
            cout << endl;
        }
        
    }
    cout << endl;
}


// write an .nde file
void write_graph(char* inname){
    FILE* input_graph = fopen(inname, "r");
    ofstream output_graph;
    output_graph.open(outname);

    // file contains number of nodes, number of edges at fist line
    // and then one edge per line
    fscanf(input_graph, "%d", &N);
    output_graph << newN << endl;

    int u, v;
    // we need to skip the first N rows in input, while we add the correct degree for the newN nodes that are not removed
    for(int i = 0; i < N; i++){
        fscanf(input_graph, "%d %d", &u, &v);
        
        if( i < MAX_N && new_names[i]!=-1){
            int curr_degree = 0;
            for(auto j : G[i]){
                if(new_names[j]!= -1)
                    curr_degree++;
            }
            output_graph << new_names[i] << " " << curr_degree << endl;
        }
    }
    
    
    while(fscanf(input_graph,"%d %d",&u, &v) != EOF){
        // check if both u,v are non-deleted, if so, add their pair to the file
        if(u<MAX_N && v<MAX_N && new_names[u]!= -1 && new_names[v] != -1){
            output_graph << new_names[u] << " " << new_names[v] << endl;
            newM++;
        }
    }

    // cout << "Output graph has " << newN << " nodes and " << newM << " edges. "<< endl;

    fclose(input_graph);
    output_graph.close();
    return;
}

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
    bool good_for_current_BCC=false;    

    // Go through all non-deleted neighbors of u
    for (auto v : G[u]) { // these are only the ones smaller than MAX_N
        // cout <<v << " is a neighbor of " << u << endl;
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (!visited[v]) {
            good_for_current_BCC = false;
            parent[v] = u;
            children++;
            find_artpts(s, v);

            // if the root finishes a child, then it basically closes an articulation point
            // if it has seen s, this is the correct BCC
            if(parent[u] == -1 && found_s){
                if(DEBUG) cout << "Root finished child " << v << " having found s" << endl;
                // it could be the second BCC for t so we need to unstack until v and the rest are bad = to_remove
                int x = cat_stack.back();
                while(x != v){ // If we are in the "right" BCC add it to the vector   
                    cat_stack.pop_back();
                    x = cat_stack.back();
                }
                cat_stack.pop_back(); 

                // the rest of the nodes in the stack are guaranteed to be 'bad' = to be removed, EXCEPT FOR TARGET
                while(!cat_stack.empty()){
                    x = cat_stack.back();
                    to_remove[x]=true;
                    cat_stack.pop_back();
                }

                to_remove[t] = false;
                
                // WE SHOULD BE ABLE TO RETURN HERE!
                return;
            }

            // if root finishes a child without finding s, all nodes in the stack are to be removed 
            if(parent[u]==-1 && !found_s){
                if(DEBUG) cout << "Root finished child "<< v <<" without finding s" << endl;
                int x = cat_stack.back();
                while(x != v){  
                    cat_stack.pop_back();
                    to_remove[x]=true;
                    if(DEBUG) cout << "Marking node " << x << " to be removed" << endl;
                    x = cat_stack.back();
                }
                cat_stack.pop_back(); 
                to_remove[v]=true;
                if(DEBUG) cout << "Marking node " << v << " to be removed" << endl;
            }

            // Check if the subtree rooted with v has
            // a connection to one of the ancestors of u
            low[u] = min(low[u], low[v]);

            // If u is not root and low value of one of
            // its child is more than discovery value of u.
            if (parent[u] != -1 && low[v] >= disc[u]){ // here is where I close an articulation point
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
                    // remove stuff from stack without adding to to_remove
                    int x = cat_stack.back();
                    while(x != v){ // If we are in the "right" BCC add it to the vector                            
                        cat_stack.pop_back();
                        x = cat_stack.back();
                    }
                    cat_stack.pop_back(); 

                    current_s = u; 
                }
                    
                if(!good_for_current_BCC){ // if I am not good, I need to delete all nodes in the stack with discovery time >=v
                    // remove stuff and delete
                    int x = cat_stack.back();
                    while(x != v){ 
                        cat_stack.pop_back();
                        to_remove[x] = true;
                        x = cat_stack.back();
                    }
                    cat_stack.pop_back();
                    to_remove[v]=true;
                }
            }
        }

        // Update low value of u for parent function calls.
        else if (v != parent[u])
            low[u] = min(low[u], disc[v]);
    
    }
    
    return;
}

// start visit for finding articulation points from t
void find_caterpillar(int s, int t)
{   
    cat_stack.erase(cat_stack.begin(), cat_stack.end());
    visit_time = 0;
    found_s = false;
    current_s = s;
    last_art = -1;


    for(int i = 0; i < G.size(); i++){
        visited[i] = false;
        parent[i] = -2;
        to_remove[i] = false;
        new_names[i] = -1;
    }
    
    parent[t] = -1;

    // cout << "About to find art points from "<< s << " to "<< t <<endl;
    // if(DEBUG) printGraph();
    // uint64_t start = timeMs();

     // only interested in the ones from s to t = caterpillar
    find_artpts(s, t);

    int curr_index = 0;
    for(int i = 0; i<G.size(); i++){
        if(visited[i] && !to_remove[i]){ // if the node is to be kept 
            new_names[i] =  curr_index++;
            newN++; // also compute the new number of nodes
            // cout << "Keeping " << i << endl;
        }
    }

    if(DEBUG){
        cout << "New names: "<< endl;
        for(int i =0; i<new_names.size(); i++){
            if(new_names[i]!=-1){
                cout << i << " -> " << new_names[i] << endl;
            }
        }

        cout << endl;
        // for each node, look through its adjacency list
        for(int i =0; i < G.size(); i++)
        {
            if(new_names[i]!=-1){
                cout << i << ": ";
                for (int j = 0; j< G[i].size(); j++){
                    if(new_names[G[i][j]]!= -1)
                        cout << '(' << i << ", " << G[i][j] << "); ";
                }
                cout << endl;
            }
            
        }
        cout << endl;
    }

    return;
}


// recursive DFS procedure from node u
void DFS(int u){
    reachable[u] = 1;

    for(int i = 0; i < G[u].size(); i++)
    {
        // reachable being equal to zero means that the node has not been deleted nor already visited
        if(!reachable[G[u][i]])
            DFS(G[u][i]);
    }

    return;
}


// iterative DFS procedure from node u
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


int main(int argc, char** argv) { 
    if(argc < 3){
        cout << "USAGE: " << argv[0] << " <input filename> <MAX_N> [output-directory]\n";
        return 0;
    }

    MAX_N = atoi(argv[2]);
    
    create_graph(argv[1]);
    
    // initialize all nodes as non-reachable
    for(int i = 0; i< reachable.size(); i++)
        reachable[i] = 0;

    srand (time(NULL));
    // srand(42);

    s = rand() % G.size();
    t = rand() % G.size();

    if(DEBUG) cout << "Checking s=" << s << " and t=" << t << endl;

    // chech reachability of s from t
    if(s!= t)
        DFS(t);

    while(!reachable[s]){
        s =  rand() % G.size();
        t = rand() % G.size();
        // cout << "Checking s=" << s << " and t=" << t << endl; 
        if(s!= t){
            for(int i = 0; i< reachable.size(); i++)
                reachable[i] = 0;

            DFS(t);
        }
    }
    
    find_caterpillar(s, t);

    if(!visited[s])
        throw logic_error("s is not reachable from t even after DFS check!");


    if(argc >= 4) outname = argv[3];
    else outname = "./"; // current directory if not specified
    
    string inname = argv[1];

    // int inlen = sizeof(argv[1])/sizeof(*argv[1]);
    // cout << "size is " << inname.length() << endl;

    // outname = "./preprocessed_datasets/";
    // cout << "outname=" << outname << endl << flush;

    for(int i=7; i<inname.length()-4; i++)
        outname = outname+(inname[i]);

    outname = outname + "-"+ to_string(new_names[s])+"-"+to_string(new_names[t])+".nde";
    
    if(DEBUG) cout << "Outname is " << outname << endl<<flush;

    
    write_graph(argv[1]);

    cout << outname << " " << new_names[s] << " " << new_names[t] << endl;

    return 0;
}
