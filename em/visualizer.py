import numpy as np
from glumpy import app, gloo, gl


import molecule
import reader

class Visualizer():

    vertex = """
        uniform vec2 resolution;
        attribute vec3 center;
        attribute float radius;
        varying vec3 v_center;
        varying float v_radius;
        void main()
        {
            v_radius = radius;
            v_center = center;
            gl_PointSize = 2.0 + ceil(2.0*radius);
            gl_Position = vec4(2.0*center.xy/resolution-1.0, v_center.z, 1.0);
        } """

    fragment = """
        varying vec3 v_center;
        varying float v_radius;
        void main()
        {
            vec2 p = (gl_FragCoord.xy - v_center.xy)/v_radius;
            float z = 1.0 - length(p);
            if (z < 0.0) discard;

            gl_FragDepth = 0.5*v_center.z + 0.5*(1.0 - z);

            vec3 color = vec3(0.7, 0.7, 0.7);
            vec3 normal = normalize(vec3(p.xy, z));
            vec3 direction = normalize(vec3(1.0, 1.0, 1.0));
            float diffuse = max(0.0, dot(direction, normal));
            float specular = pow(diffuse, 24.0);
            gl_FragColor = vec4(max(diffuse*color, specular*vec3(1.0)), 1.0);
        } """

    def show(self, molecule):
        #Get molecule data
        z, y, x = molecule.shape()
        axis_x = np.linspace(0,x-1,x)
        axis_y = np.linspace(0,y-1,y)
        axis_z = np.linspace(0,1,z)

        Z = molecule.data().reshape(-1)
        max_z = np.max(Z)
        min_z = np.min(Z)
        Z = (Z - min_z)/(max_z-min_z)
        mean_z =np.mean(Z)
        Z[Z<mean_z]=0
        Z= Z*25

        indices = np.empty((420,3),dtype=int) 
        indices[...,0] = np.arange(x)
        indices[...,1] = np.arange(y)
        indices[...,2] = np.arange(z)

        indices = np.stack(np.meshgrid(axis_x, axis_y, axis_z), -1).reshape(-1, 3)

        V = np.zeros(x*y*z, [("center", np.float32, 3),
                   ("radius", np.float32, 1)])
        V["center"] = indices
        V["radius"] = Z.reshape(-1)

        window = app.Window(512, 512, color=(1,1,1,1))
        points = gloo.Program(self.vertex, self.fragment)
        points.bind(V.view(gloo.VertexBuffer))

        @window.event
        def on_resize(width, height):
            points["resolution"] = width, height
            
        @window.event
        def on_draw(dt):
            window.clear()
            gl.glEnable(gl.GL_DEPTH_TEST)
            points.draw(gl.GL_POINTS)
            
        app.run()

  

        
filename = "../../emd_2847.map"
myreader = reader.Reader(filename)
myMolecule = myreader.read()
v = Visualizer()
v.show(myMolecule)