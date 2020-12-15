# source https://www.cs.mcgill.ca/~aassaf9/python/algorithm_x.html
import numpy as np


def solve(X, Y, count_update=[], solution=[]):
  if not X:
    yield list(solution)
    if count_update != []:
      count_update[0] += 1
  else:
    c = min(X, key=lambda c: len(X[c]))
    for r in list(X[c]):
      solution.append(r)
      ### Count updates
      if count_update != []:
        count_update[0] += 1
      cols = select(X, Y, r, count_update=count_update)
      for s in solve(X, Y, count_update=count_update, solution=solution):
        yield s
        if count_update != []:
          count_update[0] += 1
      deselect(X, Y, r, cols, count_update=count_update)
      solution.pop()
      ### Count updates
      if count_update != []:
        count_update[0] += 1


def select(X, Y, r, count_update=[]):
  cols = []
  for j in Y[r]:
    for i in X[j]:
      for k in Y[i]:
        if k != j:
          X[k].remove(i)
          ### Count updates
          if count_update != []:
            count_update[0] += 1
    cols.append(X.pop(j))
    ### Count updates
    if count_update != []:
      count_update[0] += 1
  return cols


def deselect(X, Y, r, cols, count_update=[]):
  for j in reversed(Y[r]):
    X[j] = cols.pop()
    ### Count updates
    if count_update != []:
      count_update[0] += 1
    for i in X[j]:
      for k in Y[i]:
        if k != j:
          X[k].add(i)
          ### Count updates
          if count_update != []:
            count_update[0] += 1


def gen_x_dicc(Y, initial_matrix):
  X = {i for i in range(len(initial_matrix[0]))}

  X = {j: set() for j in X}
  for i in Y:
    for j in Y[i]:
      X[j].add(i)

  return X


def gen_y_dicc(initial_matrix):
  Y = {}
  for i in range(len(initial_matrix)):
    Y[i] = np.where(np.array(initial_matrix[i]) == 1)[0].tolist()

  return Y
