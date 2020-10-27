'''
Created on Jul 19, 2020

@author: luis98
'''
from general_utils.string_utils import get_float_between_ss, get_float_between, get_vector_between, get_matrix_between


class FitMapResult(object):
  '''
  classdocs
  '''

  def __init__(self, text, type_c):
    # print(text)
    if type_c == 'n':
      self.num_points = int(get_float_between(text, "using", " points"))
      self.correlation = get_float_between(text, "correlation =", ", correlation about")
      self.correlation_about_mean = get_float_between(text, ", correlation about mean =", ", overlap")
      self.overlap = get_float_between(text, " overlap =", "\n  steps =")
      self.steps = int(get_float_between(text, "steps =", ", shift = "))
      self.shift = get_float_between(text, " shift =", ", angle = ")
      self.angle = get_float_between(text, ", angle =", " degrees\n\n")
      self.matrix_rt = get_matrix_between(text, "  Matrix rotation and translation\n", "\n  Axis")
      self.axis = get_vector_between(text, "  Axis", "\n  Axis point")
      self.axis_point = get_vector_between(text, "  Axis point", "\n  Rotation angle (degrees)")
      self.rotation_angle = get_float_between(text, "  Rotation angle (degrees)", "\n  Shift along axis")
      self.shift_along_axis = get_float_between(text, "Shift along axis", "\n\n, correlation")
    elif type_c == 'x':
      self.num_points = int(get_float_between_ss(text, "using", " points"))
      self.correlation = get_float_between(text, "correlation =", ", correlation about")
      self.correlation_about_mean = get_float_between(text, ", correlation about mean =", ", overlap")
      self.overlap = get_float_between(text, " overlap =", "\n  steps =")
      self.steps = int(get_float_between(text, "steps =", ", shift = "))
      self.shift = get_float_between(text, " shift =", ", angle = ")
      self.angle = get_float_between_ss(text, ", angle =", " degrees\n")
      '''Information about this matrix'''

      self.matrix_rt = get_matrix_between(text, "  Matrix rotation and translation\n", "\n  Axis")
      self.axis = get_vector_between(text, "  Axis", "\n  Axis point")
      self.axis_point = get_vector_between(text, "  Axis point", "\n  Rotation angle (degrees)")
      self.rotation_angle = get_float_between(text, "  Rotation angle (degrees)", "\n  Shift along axis")
      self.shift_along_axis = get_float_between_ss(text, "Shift along axis", "\n")
      self.overlap_mass = get_float_between_ss(text, "Enclosed volume for surface (#4.1) =", "\n")

    self.shape_cube1 = None
    self.shape_cube2 = None
    self.center_point1_o = None
    self.center_point2_o = None
    self.center_point1 = None
    self.center_point2 = None
    self.real_shape_cube1 = None
    self.real_shape_cube2 = None
    self.center_point1_a = None
    self.center_point2_a = None
    self.move_vector_map1 = [0, 0, 0]
    self.move_vector_map2 = [0, 0, 0]
    self.percentage_of_overlap = 0

  def get_num_points(self):
    return self.num_points

  def get_correlation(self):
    return self.correlation

  def get_correlation_about_mean(self):
    return self.correlation_about_mean

  def get_overlap(self):
    return self.overlap

  def get_steps(self):
    return self.steps

  def get_shift(self):
    return self.shift

  def get_angle(self):
    # In degrees
    return self.angle

  def get_matrix_rt(self):
    return self.matrix_rt

  def get_axis(self):
    return self.axis

  def get_axis_point(self):
    return self.axis_point

  def get_rotation_angle(self):
    # In degrees
    return self.rotation_angle

  def get_shift_along_axis(self):
    return self.shift_along_axis

  def print_data(self):
    print("Num_poins: ", self.num_points)
    print("Correlation: ", self.correlation)
    print("Correlation_about_mean: ", self.correlation_about_mean)
    print("Overlap: ", self.overlap)
    print("Steps: ", self.steps)
    print("Shift: ", self.shift)
    print("Angle: ", self.angle)
    print("Matrix_rt: ", self.matrix_rt)
    print("Axis: ", self.axis)
    print("Axis_point: ", self.axis_point)
    print("Rotation_angle: ", self.rotation_angle)
    print("Shift_along_axis: ", self.shift_along_axis)
    print("Percentage of overlap: ", self.percentage_of_overlap)
