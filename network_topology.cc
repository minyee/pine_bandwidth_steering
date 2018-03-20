#include "network_topology.h"
#include <iostream>
#include <fstream>
#include <queue>
/******************************************************************************
************************** Network node section *******************************
*******************************************************************************/

network_node::network_node(uint16_t id) : id_(id) {
};

network_node::~network_node() {
};

void network_node::add_neighbor(network_node* neighbor) {
	// neighbor is null just don't do anything and return to the caller function
	// or if the neighbor has already been added before
	if (neighbor == nullptr || (neighbors_.find(neighbor->get_id()) != neighbors_.end()))
		return;
	neighbors_.insert({neighbor->get_id(), neighbor});
};

uint16_t network_node::get_id() const {
	return id_;
};

network_node* network_node::get_neighbor(uint16_t id) const {
	network_node* target = nullptr;
	auto it = neighbors_.find(id);
	if (it != neighbors_.end()) {
		target = it->second;
	}
	return target;
};

uint16_t network_node::num_outgoing_edges() const {
	return neighbors_.size();
}

void network_node::get_all_neighbors(std::vector<network_node*>& neighbor_vector) const {
	auto it = neighbors_.begin();
	neighbor_vector.resize(0);
	while (it != neighbors_.end()) {
		neighbor_vector.push_back(it->second);
		it++;
	}
	return;
} 
/******************************************************************************
*******************************************************************************
*******************************************************************************/



/******************************************************************************
************************ Network topology section *****************************
*******************************************************************************/
network_topology::network_topology(std::string fname) : topol_filename_(fname){
	read_topol_file();
};

/**
 * Note that the topology file should take in an adjacency list, not a matrix
 **/
void network_topology::read_topol_file() {
	std::ifstream infile(topol_filename_);
	std::string line;
	
	// read in the size of the network topology first
	uint16_t size;
	infile >> size;
	nnodes_ = size;
	// now initiate all of the network nodes first
	for (int id = 0; id < size; id++) {
		nodes_.insert({id, new network_node(id)});
	}
	// now time to read the adjacency list file and then
	const char delim = ' ';
	while (std::getline(infile, line)) {
		char* dup = strdup(line.c_str());
		char* token = std::strtok(dup, &delim);
		uint16_t curr_node_id = std::atoi(token);
		auto curr_node_it = nodes_.find(curr_node_id);
		network_node* curr_node = curr_node_it->second;
		do {
		 	token = std::strtok(NULL, &delim);
		 	uint16_t neighbor_id = std::atoi(token);
		 	auto neighbor_it = nodes_.find(neighbor_id);
		 	network_node* neighbor_node = neighbor_it->second;
		 	curr_node->add_neighbor(neighbor_node);
		 } while (token != nullptr); 
	}
	infile.close();
};

network_topology::~network_topology() {
	auto it = nodes_.begin();
	while (it != nodes_.end()) {
		delete it->second;
		it++;
	}
	return;
};

network_node* network_topology::get_node(uint16_t id) const {
	network_node* found = nullptr;
	auto it = nodes_.find(id);
	if (it != nodes_.end()) {
		found = it->second;
	}
	return found;
};

void network_topology::add_node(uint16_t id) {
	auto it = nodes_.find(id);
	// the node already exists
	if (it != nodes_.end()) {
		return;
	}
	network_node* new_node = new network_node(id);
	nodes_.insert({id, new_node});
	return;
};

bool network_topology::get_outgoing_edges(uint16_t id, uint16_t& nedges) const {
	auto it = nodes_.find(id);
	if (it == nodes_.end()) {
		return false;
	}
	network_node* node = it->second;
	nedges = node->num_outgoing_edges();
	return true;
}

void network_topology::add_node_neighbor(uint16_t src, uint16_t neighbor) {
	auto src_it = nodes_.find(src);
	if (src_it == nodes_.end()) {
		add_node(src);
	}

	auto neighbor_it = nodes_.find(neighbor);
	if (neighbor_it == nodes_.end()) {
		add_node(neighbor);
	}

	auto src_node_it = nodes_.find(src);
	auto neighbor_node_it = nodes_.find(neighbor);

	(src_node_it->second)->add_neighbor(neighbor_node_it->second);
};

// perform an bfs on the topology starting from src node and then count the number of paths 
// that can potentially route src to dst.
// depth is the maximum distance DFS will traverse to account for a node, should be 2 for Flexfly-like networks
uint32_t network_topology::num_paths(uint16_t src, uint16_t dst, int depth) const {
	// if either src or dst do not exist, then just return 0 paths
	auto src_it = nodes_.find(src);
	auto dst_it = nodes_.find(dst);
	if (src_it == nodes_.end() || dst_it == nodes_.end()) {
		return 0;
	}
	std::vector<uint32_t> paths_found(nnodes_);
	std::vector<bool> visited(nnodes_); // initiate the visited nodes
	std::vector<uint8_t> distance(nnodes_);
	std::fill(paths_found.begin(), paths_found.end(), 0);
	std::fill(visited.begin(), visited.end(), false);
	std::fill(distance.begin(), distance.end(), depth * 2);
	std::queue<network_node*> q;
	q.push(src_it->second);
	network_node* curr_node = src_it->second; 
	network_node* dst_node = dst_it->second;
	distance[curr_node->get_id()] = 0;
	paths_found[curr_node->get_id()] = 1;
	while (!q.empty()) {
		curr_node = q.front();
		q.pop();
		visited[curr_node->get_id()] = true;
		std::vector<network_node *> neighbors;
		curr_node->get_all_neighbors(neighbors);
		for (auto neighbor : neighbors) {
			paths_found[neighbor->get_id()] += paths_found[curr_node->get_id()];
			distance[neighbor->get_id()] = (distance[curr_node->get_id()] + 1) < distance[neighbor->get_id()] ? (distance[curr_node->get_id()] + 1) : distance[neighbor->get_id()];
			if (!visited[neighbor->get_id()] && (distance[neighbor->get_id()] < depth)) {
				q.push(neighbor);
			}
		}
	}
	return paths_found[dst];
}

/******************************************************************************
*******************************************************************************
*******************************************************************************/