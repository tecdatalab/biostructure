from process_mrc.miscellaneou import get_mrc_segments
from process_graph.process_graph_utils import generate_graph
from process_graph.graph_algorithm import graph_aligning
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib

segments = get_mrc_segments("../../maps/1010/EMD-1010.map", 7, 3, 1)
graph = generate_graph(segments, 3, 0, 6, 1) #Preguntar con los no conectados y sub grafos
result = graph_aligning(graph, graph, 0.001)
print(result)

matplotlib.use('TKAgg')
nx.draw(graph, with_labels=True)
plt.show()