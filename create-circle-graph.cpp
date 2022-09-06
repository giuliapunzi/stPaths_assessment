#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;
string outname;
int N,M;
vector<int> edge_added;


// script to create a file containing a random sparse graph of size (both nodes and edges) given in input
void create_circle(){
    ofstream output_file; 
    output_file.open(outname);
    output_file << N << endl;
    for(int i=0; i<N; i++)
        output_file << 1 << " " << 1 << endl; // dummy rows for nde

    // add circle
    for (int i = 0; i < N-1; i++){
        output_file << i << " " << i+1 << endl;
    }
    output_file << N-1 << " " << 0 << endl;
    
    
    // for M times, extract a random pair in (0,N-1) of non-neighboring, non-covered nodes
    for(int i = 0; i < M ; i++){
        int u = rand() % N; // u in the range 0 to N-1
        int v = rand() % N;
        vector<int>::iterator find_u = find(edge_added.begin(), edge_added.end(),u);
        while(find_u !=edge_added.end()) {
            u = rand() % N;
            find_u = find(edge_added.begin(), edge_added.end(),u);
        }

        // cout << "Chosen u is " << u<< endl;

        vector<int>::iterator find_v = find(edge_added.begin(), edge_added.end(),v);
        while(v == u || abs(v-u) ==1 || abs(v-u) == N-1 || find_v != edge_added.end()){ // find  v that is not neighboring nor equal to u, and that has not been covered yet
            // cout << "Chosen v was " << v << "; rerolling" << endl;
            v = rand() % N;
            find_v = find(edge_added.begin(), edge_added.end(),v);
        }

        output_file << u << " " << v << endl;
        edge_added.push_back(u);
        edge_added.push_back(v);
        
    }

    output_file.close();
}


int main(int argc, char** argv) { 
    if(argc < 3){
        cout << "USAGE: " << argv[0] << " <num-nodes> <num-edges> [TARGET-FILENAME]\n";
        return 0;
    }


    N = atoi(argv[1]);
    M = atoi(argv[2]);

    if(M>N/2)
        throw invalid_argument("In 3-regular circle graph number of edges must be at most half the number of nodes!");

    if(argc >= 4) outname = argv[3];
    else outname = "circle-"+ to_string(N)+"-"+to_string(M)+".nde";
    

    srand (time(NULL));
    create_circle();

    return 0;
}