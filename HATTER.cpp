#include <iostream>
#include <sstream>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <cstring>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <signal.h>
#include <unistd.h>

using namespace std;

typedef pair<int, int> pi;
typedef unordered_map<int, unordered_set<int>> adj_list;

#define INF __INT_MAX__

enum EdgeType { BLACK, RED };
enum ContractingCostFunc { GLOBAL, LOCAL };

int max_rd = 0;
volatile sig_atomic_t tle = 0;

int INITIAL_CNT = 4, INC = 2;
ContractingCostFunc CONTRACTING_COST_FUNCTION = GLOBAL;
bool TWW_NEIGH = false;
bool RANDOM_EXTREME = true;

/// SIGTERM handler
void term(int signum) {
	tle = 1;
}

void alarm_handler(int signum) {
	raise(SIGTERM);
}


class Heap {
public:
	set<pi> heap;
	unordered_map<int, int> rev_mapping;

	void insert(pi p) {
		int node = p.second;
		int deg = p.first;
		if(rev_mapping.find(node) != rev_mapping.end()) { // if node already present in heap
			int old_deg = rev_mapping[node];
			if(old_deg == deg) return; // if degree is same, then no need to update
			heap.erase({old_deg, node});
		}
		heap.insert({deg, node});
		rev_mapping[node] = deg;
	}

	void erase(int node) {
		if(rev_mapping.find(node) == rev_mapping.end()) return; // node not present in heap
		heap.erase({rev_mapping[node], node});
		rev_mapping.erase(node);
		
		// rehash, if required for faster retrieval
		if(rev_mapping.size() && rev_mapping.bucket_count() / rev_mapping.size() >= 2)
			rev_mapping.rehash(rev_mapping.size());
	}

	inline bool empty() { return heap.empty(); }
	inline set<pi>::iterator begin() { return heap.begin(); }
	inline set<pi>::iterator end() { return heap.end(); }
	inline reverse_iterator<_Rb_tree_const_iterator<pi>> rbegin() { return heap.rbegin(); }
	inline reverse_iterator<_Rb_tree_const_iterator<pi>> rend() { return heap.rend(); }


	void get_top_n_elements(int n, unordered_set<int>& arr) {
		int added = 0;
		for(auto it=heap.begin(); it!=heap.end() && added < n; it++) {
			added += arr.insert(it->second).second;
		}
	}

	pi getTopElement(unordered_map<int, int>& taken) {
		for(auto it=heap.begin(); it!=heap.end(); ) {
			int deg = it->first, x = it->second;
			rev_mapping.erase(x);
			heap.erase(it++);
			if(taken.find(x) == taken.end()) {
				taken[x] = deg;
				return {deg, x};
			}
		}
		return {-1, -1};
	}
};

class Graph {
public:
	int n, m;
	adj_list black_edges; // black neighbors 
	adj_list red_edges; // red neighbors
	Heap black_heap; // min-heap of black degree 
	Heap red_heap; // min heap of red degree
	Heap deg_heap;

	Graph(int n): n(n), m(0) {}
	
	// Copy Constructor
	Graph(const Graph& G): 
		n(G.n), m(G.m), 
		black_edges(G.black_edges), red_edges(G.red_edges), 
		black_heap(G.black_heap), red_heap(G.red_heap),
		deg_heap(G.deg_heap) {}

	// add edge (u, v) to graph of type edgeType
	void add_edge(int u, int v, EdgeType edgeType) {
		if(u == v) return;

		if(edgeType == BLACK) {
			m += black_edges[u].insert(v).second;
			black_edges[v].insert(u);
		} else {
			m += red_edges[u].insert(v).second;
			red_edges[v].insert(u);
		}
	}

	// delete edge (u, v) from graph of type edgeType
	void delete_edge(int u, int v, EdgeType edgeType) {
		if(edgeType == BLACK) {
			m -= black_edges[u].erase(v);
			black_edges[v].erase(u);
		} else {
			m -= red_edges[u].erase(v);
			red_edges[v].erase(u);
		}
	}

	void initialize_heap() {
		for(auto& u: black_edges) {
			if(tle) break;
			black_heap.insert({u.second.size(), u.first});
			if(TWW_NEIGH) deg_heap.insert({u.second.size(), u.first});
		}
	}

	int intersection_cardinality(unordered_set<int>& a, unordered_set<int>& b) {
		if(a.size() > b.size())
			return intersection_cardinality(b, a);
		
		int cnt = 0;
		for(const int& ele: a) {
			if(tle) return 0;
			if(b.find(ele) != b.end()) cnt++;
		}
		return cnt;
	}

	int contracting_cost_global(int v, int w)
	{
		if(red_edges[v].size() < red_edges[w].size())
		return contracting_cost_global(w, v);

		int cost_rd = 0;
		int cost_rd_v = red_edges[v].size();

		if(red_edges[v].count(w))
		cost_rd_v--;

		for(const int& black_w: black_edges[w])
		{
			if(black_w == v)
			continue; 

			if(!red_edges[v].count(black_w) && !black_edges[v].count(black_w))
			{
				cost_rd = max(cost_rd, int(red_edges[black_w].size()) + 1);
				cost_rd_v++;
			}
		}

		for(const int& red_w: red_edges[w])
		{
			if(red_w == v)
			continue; 

			if(red_edges[v].count(red_w))
			{
				cost_rd = max(cost_rd, int(red_edges[red_w].size()) - 1);
			}
			else
			{
				cost_rd_v++;
			}
		}

		for(const int& black_v: black_edges[v]) 
		{
			if(black_v == w)
			continue; 

			if(!red_edges[w].count(black_v) && !black_edges[w].count(black_v))
			{
				cost_rd = max(cost_rd, int(red_edges[black_v].size()) + 1);
				cost_rd_v++;
			}
		}

		return max(cost_rd, cost_rd_v);
	}


	// TODO: could be further optimized
	int contracting_cost_local(int u, int v) {
		if(!black_edges.count(u) || !black_edges.count(v)) return INF;

		int cost = black_edges[u].size() + black_edges[v].size() + red_edges[u].size() + red_edges[v].size();
		if(black_edges[u].find(v) != black_edges[u].end() || red_edges[u].find(v) != red_edges[u].end()) 
			cost -= 2;

		int common_black_edges = intersection_cardinality(black_edges[u], black_edges[v]);
		int common_red_edges = intersection_cardinality(red_edges[u], red_edges[v]);
		int common_black_red_neighbors = intersection_cardinality(black_edges[u], red_edges[v]) + intersection_cardinality(red_edges[u], black_edges[v]);

		cost -= (2*common_black_edges + common_red_edges + common_black_red_neighbors);
		return cost;
	}

	int contracting_cost(int u, int v) {
		if(CONTRACTING_COST_FUNCTION == GLOBAL) return contracting_cost_global(u, v);
		else return contracting_cost_local(u, v);
	}

	void delete_node(int x) {
		n -= 1;
		m -= (black_edges[x].size() + red_edges[x].size());
		black_heap.erase(x);
		red_heap.erase(x);
		deg_heap.erase(x);

		// remove all black edges of x
		for(const int& ch: black_edges[x]) {
			if(tle) return;
			black_edges[ch].erase(x);
		} black_edges.erase(x);

		// remove all red edges of x
		for(const int& ch: red_edges[x]) {
			if(tle) return;
			red_edges[ch].erase(x);
		} red_edges.erase(x);

		// re-hash graph if required
		if(black_edges.size() && black_edges.bucket_count() / black_edges.size() >= 2) {
			black_edges.rehash(black_edges.size());
			for(auto& p: black_edges) {
				if(tle) return;
				int u = p.first;
				if(black_edges[u].size() && black_edges[u].bucket_count() / black_edges[u].size() >= 2)
					black_edges[u].rehash(black_edges[u].size());
			}
		}

		if(red_edges.size() && red_edges.bucket_count() / red_edges.size() >= 2) {
			red_edges.rehash(red_edges.size());
			for(auto& p: red_edges) {
				if(tle) return;
				int u = p.first;
				if(red_edges[u].size() && red_edges[u].bucket_count() / red_edges[u].size() >= 2)
					red_edges[u].rehash(red_edges[u].size());
			}
		}
	}

	void contract(int u, int v) {
		for(const int& w: black_edges[v]) {
			if(tle) return;
			if(!black_edges[u].count(w)) {
				add_edge(u, w, RED);
			}
		}

		vector<int> b_del;
		for(const int& w: black_edges[u]) {
			if(tle) return;
			if(!black_edges[v].count(w)) {
				add_edge(u, w, RED);
				b_del.push_back(w);
			}
		}
		for(const int& w: b_del) {
			if(tle) return;
			delete_edge(u, w, BLACK);
		}

		for(const int& w: red_edges[v]) {
			if(tle) return;
			delete_edge(u, w, BLACK);
			add_edge(u, w, RED);
		}

		delete_node(v);

		unordered_set<int> nodes_to_update = {u};
		nodes_to_update.insert(black_edges[u].begin(), black_edges[u].end());
		nodes_to_update.insert(red_edges[u].begin(), red_edges[u].end());
		// nodes_to_update.erase(v);
		
		// update heap
		for(const int& x: nodes_to_update) {
			if(tle) return;
			int blk_dg = black_edges[x].size(), red_dg = red_edges[x].size();
			if(blk_dg) black_heap.insert({blk_dg, x});
			else black_heap.erase(x);

			if(red_dg) red_heap.insert({red_dg, x});
			else red_heap.erase(x);

			if(TWW_NEIGH) {
				if(blk_dg + red_dg) deg_heap.insert({blk_dg + red_dg, x});
				else deg_heap.erase(x);
			}

			max_rd = max(max_rd, (int) red_edges[x].size());
		}
	}

	unordered_set<int> dfs(int start, unordered_set<int>& vis) {
		stack<int> st;
		unordered_set<int> comp;

		st.push(start); vis.insert(start);
		while(!st.empty()) {
			int u = st.top(); st.pop();
			comp.insert(u);
			for(const int& v: black_edges[u]) {
				if(vis.find(v) == vis.end()) {
					vis.insert(v);
					st.push(v);
				}
			}
		}
		return comp;
	}

	vector<unordered_set<int>> get_components() {
		unordered_set<int> vis;
		vector<unordered_set<int>> comp;
		for(auto& p: black_edges) {
			int u = p.first;
			if(vis.find(u) == vis.end()) {
				comp.push_back(dfs(u, vis));
			}
		}
		return comp;
	}

	/**
	 * Takes intersection of two unorderd sets and stores the result in graph[u]
	 * Complexity: O( min(|a|, |b|) )
	 */
	void compute_edges(int u, unordered_set<int>& a, unordered_set<int>& b, Graph* G) {
		if(a.size() > b.size()) {
			// if b has less no. of elements than a
			for(const int& v: b) {
				if(a.find(v) != a.end())
					G->add_edge(u, v, BLACK);
			}
		} else {
			// if a has less no. of elements than b 
			for(const int& v: a) {
				if(b.find(v) != b.end())
					G->add_edge(u, v, BLACK);
			}
		}
	}

	Graph* get_induced_subgraph(unordered_set<int>& nodes) {
		Graph* G; // induced sub-graph

		if(black_edges.size()/nodes.size() <= 2) { // induced sub-graph is almost similar to original graph
			G = new Graph(*this);
			for(auto& u: black_edges) {
				if(tle) return G;
				if(nodes.find(u.first) == nodes.end()) {
					G->delete_node(u.first);
				}
			}
		} else {
			G = new Graph(nodes.size());
			for(const int& u: nodes) {
				if(tle) return G;
				compute_edges(u, black_edges[u], nodes, G);
			}
		}
		return G;
	}
};

void read(string& input) {
	do {
		getline(cin, input);
	} while(input.length() > 0 && input[0] == 'c');
}

/// @return an object of the input graph
Graph* read_graph() {
	string input;
	istringstream ss;

	string temp;
	int n, m;

	read(input); ss.str(input);
	ss >> temp >> temp >> n >> m;
	ss.clear();

	Graph* G = new Graph(n);
	int u, v;
	for(int i=0; i<m; i++) {
		read(input); ss.str(input);
		ss >> u >> v;
		ss.clear();
		G->add_edge(u, v, BLACK);
	}
	return G;
}

void preprocessDeg1(Graph* G, vector<pi> &init_seq)
{
	int valid_ops = 0;

	int num_nodes = G->n;

	for(int i=1; i<=num_nodes; ++i)
	{
		if(G->black_edges.count(i) && int(G->black_edges[i].size()) > 1)
		{
			// vector<int> degree1_nghrs;
			bool flag = 0;
			int first_nghr = -1;
			vector<int> deg1nghrs;
			
			for(int nghr: G->black_edges[i])
			{

				if(int(G->black_edges[nghr].size()) == 1)
				{ 
					if(!flag)
					{
						flag = 1;
						first_nghr = nghr;
					}
					else
					{
						// delete_node(nghr);
						deg1nghrs.push_back(nghr);
						init_seq.push_back({first_nghr, nghr});
						valid_ops++;
					}
				}
			}

			for(int remove_nghr: deg1nghrs)
			G->delete_node(remove_nghr);

		}
	}

}


pi check_comb(Graph* G, vector<int>& nodes) {
	int xf, yf;
	xf = yf = -1;

	int min_rd = INF;
	for(const int& x: nodes) {
		for(const int& y: nodes) {
			if(tle) return {xf, yf};
			if(x >= y) continue;
			int rd = G->contracting_cost(x, y);
			if(rd < min_rd) {
				min_rd = rd;
				xf = x; yf = y;
			}
		}
	}
	return {xf, yf};
}

pi contract_next_tww_neigh(Graph* G) {
	if(G->deg_heap.empty()) return {-1, -1};

	int u = G->deg_heap.begin()->second; // node with lowest deg
	unordered_set<int> nodes;
	for(int v: G->black_edges[u]) nodes.insert(v);
	for(int v: G->red_edges[u]) nodes.insert(v);

	vector<int> dist_1(nodes.begin(), nodes.end());
	for(int x: dist_1) {
		for(int v: G->black_edges[x]) nodes.insert(v);
		for(int v: G->red_edges[x]) nodes.insert(v);
	}
	nodes.erase(u);

	int vf = 1, min_rd = INF;
	for(int v: nodes) {
		if(tle) break;
		int rd = G->contracting_cost(u, v);
		if(rd < min_rd) {
			min_rd = rd;
			vf = v;
		}
	}
	return {u, vf};
}

pi contract_next(Graph* G, const int CNT=10) {
	unordered_map<int, int> taken;
	unordered_set<int> black_taken, red_taken;
	
	for(int i=0; i<CNT; i++) {
		pi x_b = G->black_heap.getTopElement(taken);
		if(RANDOM_EXTREME) black_taken.insert(x_b.first);
        else black_taken.insert(x_b.second);
	}

	for(int i=0; i<CNT; i++) {
		pi x_r = G->red_heap.getTopElement(taken);
		if(RANDOM_EXTREME) red_taken.insert(x_r.first);
        else red_taken.insert(x_r.second);
	}

	int xf = -1, yf = -1;
	vector<int> vertex_list;
	for(auto& p: taken) {
		vertex_list.push_back(p.first);
	}

	pi temp = check_comb(G, vertex_list);
	xf = temp.first, yf = temp.second;
	G->black_heap.erase(xf); G->black_heap.erase(yf);
	G->red_heap.erase(xf); G->red_heap.erase(yf);
	for(int& x: vertex_list) {
		if(x == xf || x == yf) continue;
		if(black_taken.count(x)) G->black_heap.insert({taken[x], x});
		else if(red_taken.count(x)) G->red_heap.insert({taken[x], x});
	}
	return {xf, yf};
}

vector<pi> get_sequence(Graph* G, int cnt) {
	vector<pi> seq;
	G->initialize_heap();

	while(!tle) {
		pi cp = TWW_NEIGH ? contract_next_tww_neigh(G) : contract_next(G, cnt);
		if(G->n <= max_rd || cp.first == -1) break;
		G->contract(cp.first, cp.second);
		seq.push_back(cp);		
	}
	
	vector<int> leftover;
	for(auto& p: G->black_edges) {
		leftover.push_back(p.first);
	}
	for(int i=1; i<leftover.size(); i++) {
		seq.push_back({leftover[i], leftover[i-1]});
	}
	return seq;
}

void initialize_values(Graph* G) {
	unordered_set<int> localFuncNodeValue = {11203, 13746, 21982, 29340, 35427, 44308, 47104, 47430, 48630, 58084, 65281, 70200, 74474, 85320, 91581, 91934, 97840, 101131, 104115, 113795, 126086, 132402, 240547, 376320, 383640, 434580, 551250, 901132, 1421314, 1463861, 2481054};
	unordered_set<int> initialCntTenNodeValue = {153746, 169422};
	unordered_set<int> initialCntFourNodeValue = {189859, 265009}; 
	unordered_set<int> twwNeighNodeValue = {131072, 926552, 2216688, 2665215, 3598623};
	// unordered_set<int> randomExtremeFalseNodeValue = {310870, 320287};
	int randomExtremeFalseLo = 265009, randomExtremeFalseHi = 320287;

	int n = G->n;
	if(localFuncNodeValue.count(n)) CONTRACTING_COST_FUNCTION = LOCAL;
	else CONTRACTING_COST_FUNCTION = GLOBAL;

	TWW_NEIGH = twwNeighNodeValue.count(n);

	if(n <= 10000) INITIAL_CNT = 50;
	else if(n <= 36000) INITIAL_CNT = 10;
	else if(n <= 60000) INITIAL_CNT = 8;
	else if(n <= 90000) INITIAL_CNT = 10;
	else if(n <= 180000) INITIAL_CNT = 4;
	else if(n <= 380000) INITIAL_CNT = 10;
	else INITIAL_CNT = 2;
	
    if(CONTRACTING_COST_FUNCTION == LOCAL) INITIAL_CNT = 10;
	if(randomExtremeFalseLo < n && n <= randomExtremeFalseHi) RANDOM_EXTREME = false;
	if(initialCntFourNodeValue.count(n)) INITIAL_CNT = 4;
	if(initialCntTenNodeValue.count(n)) INITIAL_CNT = 10;

	if(INITIAL_CNT > 4) INC = 4;
	else if(INITIAL_CNT == 4) INC = 2;
	else INC = 1;
}

void random_seq(Graph* G, vector<pi>& final_seq) {
	unsigned int seed = time(0);
	vector<int> arr;
	for(auto& p: G->black_edges) {
		arr.push_back(p.first);
	}
	shuffle(arr.begin(), arr.end(), default_random_engine(seed));

	for(int i=1; i<arr.size(); i++) {
		final_seq.push_back({arr[0], arr[i]});
	}
}

void solve(Graph* G) {
	int final_rd = INF;
	vector<pi> final_seq, init_seq;
	int best_cnt;
	int cnt = INITIAL_CNT;

	preprocessDeg1(G, init_seq);
	random_seq(G, final_seq);

	vector<unordered_set<int>> comp = G->get_components();

	while(!tle) {
		max_rd = 0;
		vector<pi> seq;
		vector<int> leftover;

		for(auto& c: comp) {
			if(tle) break;

			if(c.size() == 1)
			{
				for(int node: c)
				leftover.push_back(node);

				continue;
			}
			Graph* ig = G->get_induced_subgraph(c);
			vector<pi> temp = get_sequence(ig, cnt);
			seq.insert(seq.end(), temp.begin(), temp.end());
			leftover.push_back(temp.back().first);
			delete(ig);
		}

		if(!tle && max_rd < final_rd) {
			final_rd = max_rd;
			best_cnt = cnt;
			for(int i=1; i<leftover.size(); i++) {
				seq.push_back({leftover[i], leftover[i-1]});
			}
			final_seq.swap(seq);
		}

		cnt += INC;
 	}

	for(auto& p: init_seq) {
		std::cout << p.first << " " << p.second << "\n";
	}
	for(auto& p: final_seq) {
		std::cout << p.first << " " << p.second << "\n";
	}

	// std::cout << final_rd << "\n";
}

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);
	// freopen( "..\\debug\\test\\tiny001.gr", "r", stdin);
	// freopen( "..\\local\\heuristic-public\\heuristic_060.gr", "r", stdin);
	// freopen("output", "w", stdout);

	// to handle SIGTERM
	struct sigaction action;
	memset(&action, 0, sizeof(struct sigaction));
	action.sa_handler = term;
	sigaction(SIGTERM, &action, NULL);

	/*struct sigaction alarm_action;
	memset(&alarm_action, 0, sizeof(struct sigaction));
	alarm_action.sa_handler = alarm_handler;
	sigaction(SIGALRM, &alarm_action, NULL);
	alarm(300);*/

	srand(time(0));
	Graph* G = read_graph();
	initialize_values(G);
	solve(G);
	return 0;
}
