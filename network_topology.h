#include <unordered_map>
#include <string>
#ifndef NETWORK_TOPOLOGY
#define NETWORK_TOPOLOGY 

class network_node {
public:
	network_node(uint16_t id);
	
	~network_node();

	uint16_t get_id() const;

	void add_neighbor(network_node* nn);

	network_node* get_neighbor(uint16_t) const;

	uint16_t num_outgoing_edges() const;

	void get_all_neighbors(std::vector<network_node*>& neighbor_vector) const;
private:
	uint16_t id_;

	std::unordered_map<uint16_t, network_node *> neighbors_;
};


class network_topology {
public:
	network_topology(std::string filename);

	~network_topology();

	network_node* get_node(uint16_t id) const;

	void add_node(uint16_t id);

	void add_node_neighbor(uint16_t node, uint16_t neighbor);

	bool get_outgoing_edges(uint16_t id, uint16_t& nedges) const;

	uint32_t num_paths(uint16_t src, uint16_t dst, int depth) const;

private:
	void read_topol_file();

private:
	// the collection of nodes in the network topology
	std::unordered_map<uint16_t, network_node*> nodes_;

	std::string topol_filename_;

	uint32_t nnodes_;
};

#endif
