
class SpliceGraph:
    def __init__(self):
        self.vertices = {}
        self.edges = {}
    
    def add_edge(self, u, v):
        if u not in self.vertices:
            self.vertices[u] = {'cluster': ''}
            self.edges[u] = {}
        if v not in self.vertices:
            self.vertices[v] = {'cluster': ''}
            self.edges[v] = {}
        self.edges[u][v] = True
        self.edges[v][u] = True
    
    def add_junction(self, sqtl_data):
        # populate junction information
        phenotype_id = sqtl_data['phenotype_id']
        variant_id = sqtl_data['variant_id']
        tss_distance = int(sqtl_data['tss_distance'])
        maf = float(sqtl_data['maf'])
        ma_samples = int(sqtl_data['ma_samples'])
        ma_count = int(sqtl_data['ma_count'])
        pval_nominal = float(sqtl_data['pval_nominal'])
        slope = float(sqtl_data['slope'])
        slope_se = float(sqtl_data['slope_se'])

        junction_name = phenotype_id
        junction_strand = junction_name[-1]
        junction_coords = tuple(map(int, junction_name.split(':')[1:-1]))

        # add donor and acceptor site vertices
        donor_site = junction_coords[0] - tss_distance
        acceptor_site = junction_coords[-1] - tss_distance
        donor_vertex = f'{junction_name}_{donor_site}_{acceptor_site}_d'
        acceptor_vertex = f'{junction_name}_{donor_site}_{acceptor_site}_a'
        self.vertices[donor_vertex] = {'cluster': ''}
        self.vertices[acceptor_vertex] = {'cluster': ''}
        self.edges[donor_vertex] = {}
        self.edges[acceptor_vertex] = {}

        # add junction vertex
        junction_start, junction_end = junction_coords[0], junction_coords[-1]
        junction_id = f'{junction_name}'
        self.vertices[junction_id] = {'cluster': ''}
        self.edges[junction_id] = {}

        # connect donor and acceptor site vertices to junction vertex
        self.add_edge(donor_vertex, junction_id)
        self.add_edge(junction_id, acceptor_vertex)

        # add SQTL information to junction vertex
        sqtl_dict = {
            'variant_id': variant_id,
            'tss_distance': tss_distance,
            'maf': maf,
            'ma_samples': ma_samples,
            'ma_count': ma_count,
            'pval_nominal': pval_nominal,
            'slope': slope,
            'slope_se': slope_se
        }
        
        if 'sqtl_data' not in self.edges[junction_id]:
            self.edges[junction_id]['sqtl_data'] = {}
        
        if phenotype_id not in self.edges[junction_id]['sqtl_data']:
            self.edges[junction_id]['sqtl_data'][phenotype_id] = []
            
        self.edges[junction_id]['sqtl_data'][phenotype_id].append(sqtl_dict)


    # def add_junction(self, sqtl_data):
    #     # populate junction information
    #     phenotype_id = sqtl_data['phenotype_id']
    #     variant_id = sqtl_data['variant_id']
    #     tss_distance = int(sqtl_data['tss_distance'])
    #     maf = float(sqtl_data['maf'])
    #     ma_samples = int(sqtl_data['ma_samples'])
    #     ma_count = int(sqtl_data['ma_count'])
    #     pval_nominal = float(sqtl_data['pval_nominal'])
    #     slope = float(sqtl_data['slope'])
    #     slope_se = float(sqtl_data['slope_se'])

    #     junction_name = phenotype_id
    #     junction_strand = junction_name[-1]
    #     junction_coords = tuple(map(int, junction_name.split(':')[1:-1]))

    #     # add donor and acceptor site vertices
    #     donor_site = junction_coords[0] - tss_distance
    #     acceptor_site = junction_coords[-1] - tss_distance
    #     donor_vertex = f'{junction_name}_{donor_site}_{acceptor_site}_d'
    #     acceptor_vertex = f'{junction_name}_{donor_site}_{acceptor_site}_a'
    #     self.vertices[donor_vertex] = {'cluster': ''}
    #     self.vertices[acceptor_vertex] = {'cluster': ''}
    #     self.edges[donor_vertex] = {}
    #     self.edges[acceptor_vertex] = {}

    #     # add junction vertex
    #     junction_start, junction_end = junction_coords[0], junction_coords[-1]
    #     junction_id = f'{junction_name}'
    #     self.vertices[junction_id] = {'cluster': ''}
    #     self.edges[junction_id] = {}

    #     # connect donor and acceptor site vertices to junction vertex
    #     self.add_edge(donor_vertex, junction_id)
    #     self.add_edge(junction_id, acceptor_vertex)

    #     # add SQTL information to junction vertex
    #     '''
    #     This code adds a dictionary called sqtl_data as an attribute to the edge corresponding 
    #     to the given junction_id. The sqtl_data dictionary contains the sQTL information associated 
    #     with the junction, including the variant_id, tss_distance, maf, ma_samples, ma_count, 
    #     pval_nominal, slope, and slope_se. This information can later be retrieved from the sqtl_data 
    #     dictionary for any junction in the splice graph that has sQTL data associated with it.
    #     '''
    #     #TODO: The sqtl_data of only the last sqtl of that phenotype_id gets stored. It should store all the sqtls
    #     #TODO: contd. probably because only on node of phenotype_id is stored as vertex. 
    #     #TODO: contd. hence move the current 'sqtl_data' from edges[junction_id] to a different element.
        
    #     self.edges[junction_id]['sqtl_data'] = {
    #         'variant_id': variant_id,
    #         'tss_distance': tss_distance,
    #         'maf': maf,
    #         'ma_samples': ma_samples,
    #         'ma_count': ma_count,
    #         'pval_nominal': pval_nominal,
    #         'slope': slope,
    #         'slope_se': slope_se
    #     }
    #     # print(self.edges[junction_id]['sqtl_data']['variant_id'])

    def generate_splice_graph(self):
        nodes = []
        edges = []

        for v in self.vertices:
            node = {'id': v}
            if 'sqtl_data' in self.edges[v]:
                node['sqtl_data'] = self.edges[v]['sqtl_data']
            nodes.append(node)

        for u in self.edges:
            for v in self.edges[u]:
                if v > u:  # only add edges once
                    edge = {'source': u, 'target': v}
                    edges.append(edge)

        splice_graph = {'nodes': nodes, 'edges': edges}

        return splice_graph


    def view_splice_graph(splicegraph):
        G = nx.Graph()
        # add edges to networkx graph
        for v in splicegraph.vertices:
            for e in splicegraph.edges[v]:
                G.add_edge(v, e)

        # draw graph
        pos = nx.spring_layout(G)
        nx.draw(G, pos=pos, with_labels=True)
        plt.show()


    


  
