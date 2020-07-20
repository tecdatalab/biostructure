'''
Created on Jul 19, 2020

@author: luis98
'''
import numpy  

def get_float_betwen(tex,a,b):
    pos_a = tex.find(a)
    pos_b = tex.find(b)
    return float(tex[pos_a+len(a):pos_b])

def get_vector_betwen(tex,a,b):
    pos_a = tex.find(a)
    pos_b = tex.find(b)
    
    vector_tex = tex[pos_a+len(a):pos_b]
    result = []
    for j in vector_tex.split(" "):
        try:
            result.append(float(j))
        except:
            pass
    return numpy.array(result) 

def get_matriz_betwen(tex,a,b):
    pos_a = tex.find(a)
    pos_b = tex.find(b)
    
    matriz_tex = tex[pos_a+len(a):pos_b]
    result = []
    for i in matriz_tex.split("\n"):
        temp =[]
        for j in i.split(" "):
            try:
                temp.append(float(j))
            except:
                pass
        result.append(temp)
        
    return numpy.array(result) 

class FitMapResult(object):
    '''
    classdocs
    '''


    def __init__(self, text):
        #print(text)
        self.num_poins = int(get_float_betwen(text,"using "," points"))
        self.correlation = get_float_betwen(text,"correlation = ",", correlation about")
        self.correlation_about_mean = get_float_betwen(text,", correlation about mean = ",", overlap")
        self.overlap = get_float_betwen(text," overlap = ","\n  steps =")
        self.steps = int(get_float_betwen(text,"steps = ",", shift = "))
        self.shift = get_float_betwen(text," shift = ",", angle = ")
        self.angle = get_float_betwen(text,", angle = "," degrees\n\n")
        self.matrix_rt = get_matriz_betwen(text,"  Matrix rotation and translation\n","\n  Axis   ")
        self.axis = get_vector_betwen(text,"  Axis   ","\n  Axis point")
        self.axis_point = get_vector_betwen(text,"  Axis point   ","\n  Rotation angle (degrees)")
        self.rotation_angle = get_float_betwen(text,"  Rotation angle (degrees)   ","\n  Shift along axis")
        self.shift_along_axis = get_float_betwen(text,"Shift along axis   ","\n\n, correlation")
    
    def get_num_poins(self):
        return self.num_poins
    
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
        #In degress
        return self.angle
    
    def get_matrix_rt(self):
        return self.matrix_rt
    
    def get_axis(self):
        return self.axis
    
    def get_axis_point(self):
        return self.axis_point
    
    def get_rotation_angle(self):
        #In degress
        return self.rotation_angle
    
    def get_shift_along_axis(self):
        return self.shift_along_axis
    
    def print_data(self):
        print("Num_poins: ",self.num_poins)
        print("Correlation: ",self.correlation)
        print("Correlation_about_mean: ",self.correlation_about_mean)
        print("Overlap: ",self.overlap)
        print("Steps: ",self.steps)
        print("Shift: ",self.shift)
        print("Angle: ",self.angle)
        print("Matrix_rt: ",self.matrix_rt)
        print("Axis: ",self.axis)
        print("Axis_point: ",self.axis_point)
        print("Rotation_angle: ",self.rotation_angle)
        print("Shift_along_axis: ",self.shift_along_axis)
    
    
    
    
    
    
    
    
    