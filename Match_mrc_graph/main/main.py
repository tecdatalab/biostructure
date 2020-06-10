from process_mrc.miscellaneou import get_mrc_segments
from process_graph.process_graph_utils import generate_graph
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib

segments = get_mrc_segments("../../maps/1010/EMD-1010.map", 7, 3, 1)
graph = generate_graph(segments, 3, 0, 6, 1) #Preguntar con los no conectados y sub grafos

matplotlib.use('TKAgg')
nx.draw(graph, with_labels=True)
plt.show()