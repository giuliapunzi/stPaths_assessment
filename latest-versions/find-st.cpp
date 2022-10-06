#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <time.h>
#include <stack>

using namespace std;

vector<vector<int>> G;
vector<bool> reachable;

// checks if a given edge (u,v) belongs to graph G
inline bool is_edge(int u, int v)
{   
    // if(deleted[u] || deleted[v])
    //     return false;
    
    if(find(G[u].begin(), G[u].end(), v) == G[u].end())
        return false;
    else
        return true;
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
    reachable.resize(N, 1);

    // initialize reachable vector to be all 1 and degree to be zero
    for(int i = 0; i < G.size() ; i++){
        reachable[i] = 1;
        // current_degree[i] = 0;
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
            // current_degree[u]++;
            // current_degree[v]++;
        }
        
    }

    cout << "Input graph has " << N << " nodes and " << real_edges << " edges. "<< endl;


    fclose(input_graph);
    return;
}

// long long count_DFS= 0;

// recursive DFS procedure from node u
void DFS(int u){
    // count_DFS++;
    // cout << "Entering DFS for " << u << endl << flush;
    // if(count_DFS % 10000 == 0)
    //     cout << "*" << flush;

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




int main(int argc, char* argv[]){ 
     if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <FILENAME>" << std::endl;
        return 1;
    }

    char * input_filename = argv[1];
    create_graph(input_filename); // initialize 

    // initialize all nodes as non-reachable
    for(int i = 0; i< reachable.size(); i++)
        reachable[i] = 0;

    srand (42);
    int s =  rand() % G.size();
    int t = rand() % G.size();
    // cout << "Checking s=" << s << " and t=" << t << endl;

    // chech reachability of s from t
    if(s!= t)
        DFS_iter(t);

    while(!reachable[s]){
        s =  rand() % G.size();
        t = rand() % G.size();
        // cout << "Checking s=" << s << " and t=" << t << endl; 
        if(s!= t)
            DFS_iter(t);
    }

    cout << s << " " << t << " " << endl;
    return 0;
}