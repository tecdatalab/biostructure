import networkx as nx


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


def remove_node_by_name(graph, node_name):
  try:
    graph.remove_node(node_name)
  except:
    pass


def remove_edge_nodes(graph, node_name1, node_name2):
  try:
    graph.remove_edge(node_name1, node_name2)
  except:
    pass

  try:
    graph.remove_edge(node_name2, node_name1)
  except:
    pass


def graph_is_connect(graph):
  return nx.is_connected(graph)
