from collections import defaultdict
import os
import numpy as np              # For graph handling
import networkx as nx           # Visualizes graphs and runs algorithms such as MST
import matplotlib.pyplot as plt # Plots the graphs
from matplotlib.backends.backend_pdf import PdfPages

# Represent a weighted undirected graph.
class Graph:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file = os.path.join(script_dir, 'adj_matrix.csv')
    df = np.loadtxt(file, delimiter=',', dtype=int)
    """
    Since the focus is on the definitions of MBST, MMST, and MVST,
    this implementation does not focus on directed graphs.
    However, code for verifying matrix symmetry is included.
    """
    def __init__(self):
        self.nodes = len(self.df)  # Number of nodes
        self.node_list = list(range(self.nodes))
        self.symmetry = np.array_equal(self.df, self.df.T)
        
        if self.symmetry:
            # Get upper triangle indices excluding the diagonal
            row_index, col_index = np.triu_indices(self.nodes, k=1) 
            edge_dict = {}
            
            for i, j in zip(row_index, col_index):
                if self.df[i, j] > 0:
                    edge_dict[int(i), int(j)] = int(self.df[i, j])
                    
            # Sort edges by weight in non-descending order (Kruskal's sorting table in a list)
            self.sorted_edges = []
            for (u, v), w in sorted(edge_dict.items(), key=lambda x: x[1]):
                self.sorted_edges.append((u, v, w))

        else: 
            raise ValueError("Error. Please check CSV file.")

    def get_layout(self, mst_edges):
        """
        Create a fixed layout based on the MST edges.
        
        Parameters:
            mst_edges (list): List of edges in the MST.
        
        Returns:
            dict: Positions of nodes for consistent graph visualization.
        """

        # Graph used only for layout calculation
        g_layout = nx.Graph()
        g_layout.add_nodes_from(self.node_list)
        g_layout.add_weighted_edges_from(mst_edges)

        master_pos = nx.spring_layout(g_layout, seed=42)

        return master_pos
        
    
    def visualize(self, tree_edges, pos, title, show=True):
        """
        Visualize the graph given weighted edges, node positions, and title.
        
        Parameters:
            tree_edges (list): List of weighted edges to draw.
            pos (dict): Node positions for layout.
            title (str): Title of the plot.
            show (bool): Whether to display the plot immediately.
        
        Returns:
            matplotlib.figure.Figure: The generated figure object.
        """     
        g = nx.Graph()
        g.add_nodes_from(self.node_list)
        g.add_weighted_edges_from(tree_edges)

        edge_labels = nx.get_edge_attributes(g, 'weight')

        # Draw the graph
        plt.figure(figsize=(10, 7))

        # Draw all nodes from the original tree
        nx.draw_networkx_nodes(
            self.node_list,
            pos,
            node_size=700,
            node_color='lightblue',
            alpha=0.8
        )

        # Add labels
        nx.draw_networkx_labels(g, pos, font_size=12, font_family='sans-serif')
        
        # Draw edges from the specific tree
        nx.draw_networkx_edges(
            g,
            pos,
            edgelist=g.edges(),
            width=2,
            edge_color='gray'
        )

        # Add edge labels
        nx.draw_networkx_edge_labels(g, pos, edge_labels=edge_labels, font_size=10)

        plt.title(title, size=15)
        plt.axis('off')

        # Show figure if requested
        if show:
            plt.show()
        return plt.gcf() 


    def show_info(self):
        """
        Display information about the graph including the adjacency matrix.
        """
        print(
            "CSV file includes a symmetrical",
            self.nodes,
            "x",
            self.nodes,
            "adjacency matrix of a graph with nodes",
            self.node_list,
            ":\n",
            self.df
        )
        print('\n')

graph = Graph()

class DSU:
    def __init__(self, graph):
        self.parent = [i for i in range(graph.nodes)]
        self.rank = [0] * graph.nodes

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])  # Path compression
        return self.parent[x]

    def union(self, x, y):
        root_x = self.find(x)
        root_y = self.find(y)

        if root_x == root_y:  # Already connected
            return False
                                                     
        # Union by rank                              
        if self.rank[root_x] < self.rank[root_y]:
            self.parent[root_x] = root_y
        elif self.rank[root_x] > self.rank[root_y]:
            self.parent[root_y] = root_x
        else:
            self.parent[root_y] = root_x
            self.rank[root_x] += 1

        return True
        

    def num_components(self):
        """
        Count the number of disjoint sets currently in the graph.
        
        Returns:
            int: Number of connected components.
        """
        roots = {self.find(x) for x in range(len(self.parent))}
        return len(roots)

        return True

dsu = DSU(graph)

def update_optimal(new_tree, new_value, current_value, tree_list):
    """
    Compare a new tree's metric to the current best and update the list of optimal trees.
    
    Parameters:
        new_tree (list): The new tree edges.
        new_value (float): The metric value of the new tree.
        current_value (float): The current best metric value.
        tree_list (list): List of trees with the current best metric.
    
    Returns:
        tuple: Updated current value and list of trees.
    """
    if new_value < current_value:
        return new_value, [new_tree.copy()]
    elif new_value == current_value:
        tree_list.append(new_tree.copy())
        return current_value, tree_list
    else:
        return current_value, tree_list

def dedupe(tree_list):
    """
    Remove duplicate trees from the list.
    
    Parameters:
        tree_list (list): List of trees (each a list of edges).
    
    Returns:
        list: List of unique trees.
    """
    dedupe = {tuple(sorted(tree)) for tree in tree_list}
    return [list(tree) for tree in dedupe]

def STO(dsu, graph):
    """
    Perform Spanning Tree Optimization to find MST, MBST, MMST, and MVST.
    
    Parameters:
        dsu (DSU): Disjoint set union structure for the graph.
        graph (Graph): The graph object.
    
    Returns:
        tuple: MST total weight, MST edges, MBSTs, MMSTs, MVSTs, bottleneck, median, variance.
    """
    MST_sum_weight = 0
    MST_edges = []

    # Stores the value of the MBST for later use
    bottleneck = 0

    # Kruskal's main loop
    for u, v, w in graph.sorted_edges: 
        if dsu.union(u, v): 
            MST_sum_weight += w                 
            MST_edges.append((u, v, w))
            bottleneck = max(bottleneck, w)

    # MST weights list
    MST_weights = [w for  _, _, w in MST_edges]

    # MST is the first candidate for all spanning tree types
    MBSTs = [MST_edges.copy()]

    MMSTs = [MST_edges.copy()]
    median = np.median(MST_weights)

    MVSTs = [MST_edges.copy()]
    var = np.var(MST_weights)

    # Build an nx.Graph from the MST for fast pathfinding
    MST_graph = nx.Graph()
    MST_graph.add_weighted_edges_from(MST_edges)
    
    # Adding one edge outside the MST creates exactly one cycle
    # Create a set of MST edges for O(1) lookup
    MST_set = set(MST_edges)

    # Iterate only over non-MST edges
    excluded_edges = [edge for edge in graph.sorted_edges if edge not in MST_set]

    for edge_out in excluded_edges:
        u, v, w_out = edge_out

        # Find the unique path (the cycle) in the MST (O(N))
        path_nodes = nx.shortest_path(MST_graph, source=u, target=v)
        
        # Iterate over edges on this path (less than N) (O(k))
        for i in range(len(path_nodes) - 1):
            node1, node2 = path_nodes[i], path_nodes[i + 1]
            w_in = MST_graph[node1][node2]['weight']
            edge_in = (node1, node2, w_in)

            # Edge must be in the original MST to be removed
            if edge_in not in MST_set:
                continue

            # Create the new tree by replacing one edge
            tentative_tree = MST_edges.copy()
            tentative_tree.remove(edge_in)
            tentative_tree.append(edge_out)

            tentative_weights = [w for  _, _, w in tentative_tree]

            # For MBST (bottleneck)
            new_bottleneck = max(tentative_weights)
            bottleneck, MBSTs = update_optimal(
                tentative_tree, 
                new_bottleneck, 
                bottleneck, 
                MBSTs
            )
            MBSTs = dedupe(MBSTs)
            
            # For MMST (median)
            new_median = np.median(tentative_weights)
            median, MMSTs = update_optimal(tentative_tree, new_median, median, MMSTs)
            MMSTs = dedupe(MMSTs)
            
            # For MVST (variance)
            new_var = np.var(tentative_weights)
            var, MVSTs = update_optimal(tentative_tree, new_var, var, MVSTs)
            MVSTs = dedupe(MVSTs)

    return MST_sum_weight, MST_edges, MBSTs, MMSTs, MVSTs, bottleneck, median, var

def run():
    graph = Graph()

    graph.df = np.loadtxt(graph.file, delimiter=',', dtype=int)
    graph.nodes = len(graph.df)
    graph.node_list = list(range(graph.nodes))
    graph.symmetry = np.array_equal(graph.df, graph.df.T)

    if not graph.symmetry:
        raise ValueError(
            "Matrix is not symmetric. Please check the adjacency matrix."
        )
    
    dsu = DSU(graph)
    MST_sum_weight, MST_edges, MBSTs, MMSTs, MVSTs, bottleneck, median, var = STO(dsu, graph)
    master_pos = graph.get_layout(MST_edges)

    # Export information and visuals to PDF
    with PdfPages("graph_report.pdf") as pdf:

        # Original graph
        fig = plt.figure(figsize=(12, 6))
        gs = fig.add_gridspec(1, 2, width_ratios=[2, 1])

        # Left side: Graph
        ax0 = fig.add_subplot(gs[0])
        g = nx.Graph()
        g.add_weighted_edges_from(graph.sorted_edges)
        nx.draw(
            g,
            pos=master_pos,
            with_labels=True,
            node_color='lightblue',
            node_size=700,
            edge_color='gray',
            ax=ax0
        )               
        ax0.set_title("Original Graph", size=14)
        ax0.axis('off')

        # Right side: Text information
        ax1 = fig.add_subplot(gs[1])
        ax1.axis('off')
        text_str = f"Number of nodes: {graph.nodes}\nEdges:\n{graph.sorted_edges}"
        ax1.text(0, 1, text_str, va='top', ha='left', fontsize=12)
        
        pdf.savefig(fig)
        plt.close(fig)

        # MST
        fig = plt.figure(figsize=(12, 6))
        gs = fig.add_gridspec(1, 2, width_ratios=[2, 1])
        ax0 = fig.add_subplot(gs[0])
        g = nx.Graph()
        g.add_weighted_edges_from(MST_edges)
        nx.draw(
            g, 
            pos=master_pos, 
            with_labels=True, 
            node_color='lightblue', 
            node_size=700, 
            edge_color='gray', 
            ax=ax0
        )
        ax0.set_title("Minimum Spanning Tree", size=14)
        ax0.axis('off')
        ax1 = fig.add_subplot(gs[1])
        ax1.axis('off')
        text_str = f"MST total weight: {MST_sum_weight}\nEdges:\n{MST_edges}"
        ax1.text(0, 1, text_str, va='top', ha='left', fontsize=12)
        pdf.savefig(fig)
        plt.close(fig)

        # MBST
        fig = plt.figure(figsize=(12, 6))
        gs = fig.add_gridspec(1, 2, width_ratios=[2, 1])
        ax0 = fig.add_subplot(gs[0])
        g = nx.Graph()
        g.add_weighted_edges_from(MBSTs[0])
        nx.draw(
            g, 
            pos=master_pos, 
            with_labels=True, 
            node_color='lightgreen', 
            node_size=700, 
            edge_color='gray', 
            ax=ax0
        )
        ax0.set_title("Minimum Bottleneck Spanning Tree", size=14)
        ax0.axis('off')
        ax1 = fig.add_subplot(gs[1])
        ax1.axis('off')
        text_str = f"MBST count: {len(MBSTs)}\nBottleneck: {bottleneck}\nEdges:\n{MBSTs[0]}"
        ax1.text(0, 1, text_str, va='top', ha='left', fontsize=12)
        pdf.savefig(fig)
        plt.close(fig)

        # MMST
        for i, mmst in enumerate(MMSTs):
            fig = plt.figure(figsize=(12, 6))
            gs = fig.add_gridspec(1, 2, width_ratios=[2, 1])
            ax0 = fig.add_subplot(gs[0])
            g = nx.Graph()
            g.add_weighted_edges_from(mmst)
            nx.draw(
                g, 
                pos=master_pos, 
                with_labels=True, 
                node_color='lightcoral', 
                node_size=700, 
                edge_color='gray', 
                ax=ax0
            )
            ax0.set_title(f"Minimum Median Spanning Tree {i+1}", size=14)
            ax0.axis('off')
            ax1 = fig.add_subplot(gs[1])
            ax1.axis('off')
            text_str = f"Median: {median}\nEdges:\n{mmst}"
            ax1.text(0, 1, text_str, va='top', ha='left', fontsize=12)
            pdf.savefig(fig)
            plt.close(fig)

        # MVST
        fig = plt.figure(figsize=(12, 6))
        gs = fig.add_gridspec(1, 2, width_ratios=[2, 1])
        ax0 = fig.add_subplot(gs[0])
        g = nx.Graph()
        g.add_weighted_edges_from(MVSTs[0])
        nx.draw(
            g, 
            pos=master_pos, 
            with_labels=True, 
            node_color='khaki', 
            node_size=700, 
            edge_color='gray', 
            ax=ax0
        )
        ax0.set_title("Minimum Variance Spanning Tree", size=14)
        ax0.axis('off')
        ax1 = fig.add_subplot(gs[1])
        ax1.axis('off')
        text_str = f"Variance: {var}\nEdges:\n{MVSTs[0]}"
        ax1.text(0, 1, text_str, va='top', ha='left', fontsize=12)
        pdf.savefig(fig)
        plt.close(fig)

    print("All outputs (text + graphs) have been saved to 'graph_report.pdf'")
    

if __name__ == "__main__":
    run()