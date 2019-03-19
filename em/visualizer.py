import numpy as np
from glumpy import app, gloo, gl, glm
from skimage import measure


import molecule
import reader


class Visualizer():

    vertex = """
        uniform mat4   model;
        uniform mat4   view;       // View matrix
        uniform mat4   projection;   // Projection matrix
        attribute vec3 position;      // Vertex position
        
        

        void main()
        {
            gl_Position = projection * view * model * vec4(position,1.0);
        } """

    fragment = """
        void main()
        {
            gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);;
        } """

    

    def show(self, molecule):

        window = app.Window(512, 512, color=(1,1,1,1))
        phi, theta = 0, 0



        @window.event
        def on_draw(dt):
            nonlocal phi, theta
            window.clear()
            gl.glDisable(gl.GL_BLEND)
            gl.glEnable(gl.GL_DEPTH_TEST)
            gl.glEnable(gl.GL_POLYGON_OFFSET_FILL)
            
            points.draw(gl.GL_TRIANGLES, I)

            model = np.eye(4, dtype=np.float32)
            glm.rotate(model, theta, 1, 0, 0)
            glm.rotate(model, phi, 0, 1, 0)
            model = glm.translate(model, -70)
            points['model'] = model


        
        @window.event
        def on_resize(width, height):
            points['projection'] = glm.perspective(45.0, width / float(height), 1.0, 140.0)

        @window.event
        def on_mouse_drag(x, y, dx, dy, buttons):
            nonlocal phi, theta
            phi += dx
            theta += dy
        
            

        Z = molecule.data()



        verts, faces, normals, values = measure.marching_cubes_lewiner(Z)
       

        V = np.zeros(len(verts), [("position", np.float32, 3)])
        V["position"] = verts 

        V = V.view(gloo.VertexBuffer)
        I = (faces).astype(np.uint32)
        I = I.view(gloo.IndexBuffer)
               
        points = gloo.Program(self.vertex, self.fragment)
        points.bind(V)

        model = np.eye(4, dtype=np.float32)
        model = glm.translate(model, -70)
        points['model'] = model

        view = np.eye(4, dtype=np.float32)
        view = glm.translate(view, 0, 0, -140)
        points['view'] = view

        

        app.run()


filename = "tests/EMD-2677.map"
myreader = reader.Reader(filename)
myMolecule = myreader.read()
v = Visualizer()
v.show(myMolecule)