#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;
string outname;
int N,M;
vector<int> edge_added;
vector<int> degree;
vector<vector<int>> G;




// script to create a file containing a random sparse graph of size (both nodes and edges) given in input
void create_ladder(){
    ofstream output_file; 
    output_file.open(outname);
    output_file << N << endl;

    degree.resize(N);
    G.resize(N);

    for(int i=2; i < N-2; i++)
        degree[i]=3; // initialize degrees
    
    degree[0] = 2;
    degree[1]=2;
    degree[N-1] = 2;
    degree[N-2] = 2;
    
    for(int i = 0; i < N ; i+=2){
        // only add to smallest index adj list, so that no duplicates
        if(i<N-2) G[i].push_back(i+2);
        G[i].push_back(i+1);
        if(i<N-2) G[i+1].push_back(i+3);
    }

    // add degree lines
    for(int i=0; i<N; i++){
        output_file << i << " " << degree[i]<<endl;
    }    


    // add edges
    for(int u = 0; u < N; u++){
        for(auto v : G[u])
            output_file << u << " " << v << endl;
    }


    output_file.close();
    return;
}


int main(int argc, char** argv) { 
    if(argc < 2){
        cout << "USAGE: " << argv[0] << " <num-nodes> [TARGET-FILENAME]\n";
        return 0;
    }


    N = atoi(argv[1]);

    if(N % 2 != 0)
        throw invalid_argument("Number of nodes cannot be odd!");

    if(argc >= 3) outname = argv[2];
    else outname = "ladder-"+ to_string(N)+".nde";
    

    srand (time(NULL));
    create_ladder();

    return 0;
}