def get_node_by_position(graph, pos):
  all_nodes_name = list(graph.nodes)
  result = graph.nodes[all_nodes_name[pos]]
  return result


def get_node_name_by_pos(graph, pos):
  g_nodes = list(graph.nodes)
  return g_nodes[pos]


def get_node_index(graph, name):
  g_nodes = list(graph.nodes())
  return g_nodes.index(name)
