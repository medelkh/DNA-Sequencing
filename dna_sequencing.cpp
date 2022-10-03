#include <bits/stdc++.h>

#define ll long long
#define pb push_back
#define x first
#define y second
#define sz(u) (int)(u.size())

using namespace std;


int sequence_length;
int num_reads;
int read_length;
int kmer_length;    
double expected_coverage;

struct Graph{
    vector<vector<int>> adj;
    vector<string> noeuds;
    Graph(vector<vector<int>> G, vector<string> nodes): adj(G), noeuds(nodes){}
};

int gen_rand(int l, int r){
    return l + rand()%(r-l+1);
}

int alignement_lent(string a, string b){ //calculates Levenshtein distance
    int n = a.length(), m = b.length();
    vector<vector<int>> dp(n+1,vector<int>(m+1));
    dp[0][0]=0;
    for(int i=0;i<=n;i++){
        for(int j=0;j<=m;j++){
            if(i==0 && j==0) continue;
            dp[i][j] = n+m;
            int non_egaux = a[i]!=b[j];
            if(i>0) dp[i][j] = dp[i-1][j]+1;
            if(j>0) dp[i][j] = min(dp[i][j], dp[i][j-1]+1);
            if(i>0 && j>0) dp[i][j] = min(dp[i][j], dp[i-1][j-1]+non_egaux);
        }
    }
    return dp[n][m];
}

Graph  de_bruijn_graph(vector<string> lectures, int k){ //k : kmer size
    unordered_map<string,int> indice_du_kmer;
    vector<string> noeud;
    vector<vector<int>> graph_db;
    map<pair<int,int>,int> edges;
    for(string t : lectures){
        int indice_precedent = -1;
        for(int i=0;i+k<=t.length();i++){
            string kmer="";
            for(int j=0;j<k;j++){
                kmer += t[i+j];
            }
            int indice = -1;
            if(indice_du_kmer.count(kmer) == 0){
                indice = (int)(noeud.size());
                indice_du_kmer[kmer] = indice;
                graph_db.push_back(vector<int>{});
                noeud.pb(kmer);
            }
            else{
                indice = indice_du_kmer[kmer];
            }
            if(indice_precedent!=-1){
                edges[{indice_precedent,indice}]++;
                 }
            indice_precedent = indice;
        }
    }
    for(auto &p : edges){
        p.y = max(1,(int)(round(double(p.y)/expected_coverage)+0.5));
        for(int k=0;k<p.y;k++){
            graph_db[p.x.x].pb(p.x.y);
        }
    }
    return Graph(graph_db, noeud);
}


vector<int> parcours_eulerien(vector<vector<int>> graph){ //returns an eulerian path
    vector<int> parcours_temp;
    vector<int> ret;
    int init = 0, n = graph.size();
    vector<int> in_deg(n);
    for(int i=0;i<n;i++){
        for(int j=0;j<sz(graph[i]);j++){
            in_deg[graph[i][j]]++;
        }
    }
    for(int i=0;i<n;i++){
        if(sz(graph[i])==in_deg[i]+1){
            init = i;
            break;
        }
    }
    parcours_temp.pb(init);
    while(!parcours_temp.empty()){
        int cur = parcours_temp.back();
        bool sort = false;
        while(graph[cur].empty()){
            ret.pb(cur);
            parcours_temp.pop_back();
            if(parcours_temp.empty()){
                sort = true;
                break;
            }
            cur = parcours_temp.back();
        }
        if(sort) break;
        if(!graph[cur].empty()){
            parcours_temp.pb(graph[cur].back());
            graph[cur].pop_back();
        }
    }
    reverse(ret.begin(),ret.end());
    return ret;
}

void calc_cvg(){
    expected_coverage = num_reads * double(read_length - kmer_length) / double(sequence_length-read_length+1);
}
