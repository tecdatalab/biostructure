from process_mrc.miscellaneou import get_mrc_segments, get_mrc_synthetic_segments_pdb
from process_graph.process_graph_utils import generate_graph
from process_graph.graph_algorithm import graph_aligning


segments = get_mrc_segments("../../maps/1010/EMD-1010.map", 7, 3, 1)
#segments = get_mrc_synthetic_segments_pdb("../pdb_mrc/exit_pdb/175d", 7)
graph = generate_graph(segments, 50, 0, 6, 1) #Preguntar con los no conectados y sub grafos
result = graph_aligning(graph, graph, 1, False)
print(result)

#matplotlib.use('TKAgg')
#nx.draw(graph, with_labels=True)
#plt.show()