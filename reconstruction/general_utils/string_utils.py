import numpy as np


def get_float_between_ss(tex, a, b):
  pos_a = tex.find(a)
  tex = tex[pos_a + len(a):-1]
  pos_b = tex.find(b)

  return float(tex[0:pos_b])


def get_float_between(tex, a, b):
  pos_a = tex.find(a)
  pos_b = tex.find(b)
  return float(tex[pos_a + len(a):pos_b])


def get_float_value(tex, a, b):
  r = 0
  r = get_float_between_ss(tex, a, b)
  return r


def get_matrix_between(tex, a, b):
  pos_a = tex.find(a)
  pos_b = tex.find(b)

  matrix_tex = tex[pos_a + len(a):pos_b]
  result = []
  for i in matrix_tex.split("\n"):
    temp = []
    for j in i.split(" "):
      try:
        temp.append(float(j))
      except ValueError:
        pass
    result.append(temp)

  return np.array(result)


def get_vector_between(tex, a, b):
  pos_a = tex.find(a)
  pos_b = tex.find(b)

  vector_tex = tex[pos_a + len(a):pos_b]
  result = []
  for j in vector_tex.split(" "):
    try:
      result.append(float(j))
    except:
      pass
  return np.array(result)


def change_string(ip, jp, initial, final):
  result = ""
  posfinal = 0

  for i in range(len(initial)):
    if i >= ip and i <= jp:
      result += final[posfinal]
      posfinal += 1
    else:
      result += initial[i]

  return result
