import numpy as np
from glumpy import app, gloo, gl, glm
from skimage import measure


import molecule
import reader

class Visualizer():

    ''' OLD
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
            gl_PointSize = 7.0 + ceil(7.0*radius);
            gl_Position = vec4(7.0*center.xy/resolution-1.0, v_center.z, 1.0);
        } """

    fragment = """
        varying vec3 v_center;
        varying float v_radius;
        void main()
        {
            if (v_radius == 0.0) discard;

            else{
                vec2 p = (gl_FragCoord.xy - v_center.xy)/v_radius;
                vec3 color = vec3(0.7, 0.7, 0.7);
                gl_FragColor = vec4(color, 1.0);
            }
        } """

    def show(self, molecule):
        #Get molecule data
        z, y, x = molecule.shape()
        axis_x = np.linspace(0,x-1,x)
        axis_y = np.linspace(0,y-1,y)
        axis_z = np.linspace(0,z-1,z)

        Z = molecule.data().reshape(-1)
        max_z = np.max(Z)
        min_z = np.min(Z)
        Z = (Z - min_z)/(max_z-min_z)
        mean_z =np.mean(Z)
        Z[Z <mean_z] = 0
        print(np.count_nonzero(Z))

        indices = np.empty((140,3),dtype=int) 
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

  

        
filename = "tests/EMD-2677.map"
myreader = reader.Reader(filename)
myMolecule = myreader.read()
v = Visualizer()
v.show(myMolecule)
'''
########################################################

    vertex = """
        uniform mat4   model;         // Model matrix
        uniform mat4   view;          // View matrix
        uniform mat4   projection;    // Projection matrix
        attribute vec3 position;      // Vertex position
        uniform vec2 resolution;
        

        void main()
        {
            gl_Position = projection * view * model * vec4(position,1.0);
            //gl_Position = vec4(7.0*position.xy/resolution-1.0, position.z, 1.0);
        } """

    fragment = """
        void main()
        {
            gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);;
        } """

    

    def show(self, molecule):

        window = app.Window(512, 512, color=(1,1,1,1))

        @window.event
        def on_draw(dt):
            window.clear()
            gl.glDisable(gl.GL_BLEND)
            gl.glEnable(gl.GL_DEPTH_TEST)
            gl.glEnable(gl.GL_POLYGON_OFFSET_FILL)
            points['ucolor'] = .75, .75, .75, 1
            points.draw(gl.GL_TRIANGLES, I)



        @window.event
        def on_resize(width, height):
            points['projection'] = glm.perspective(45.0, width / float(height), 2.0, 100.0)
            

        Z = molecule.data()

        verts, faces, normals, values = measure.marching_cubes_lewiner(Z)



        V = np.zeros(len(verts), [("position", np.float32, 3)])
        V["position"] = verts
        V = V.view(gloo.VertexBuffer)
        I = verts[faces].astype(np.uint32)
        I = I.view(gloo.IndexBuffer)
               
        points = gloo.Program(self.vertex, self.fragment)
        points.bind(V)
        points['model'] = np.eye(4, dtype=np.float32)
        points['view'] = glm.translation(0, 0, -5)

        app.run()


filename = "tests/EMD-2677.map"
myreader = reader.Reader(filename)
myMolecule = myreader.read()
v = Visualizer()
v.show(myMolecule)