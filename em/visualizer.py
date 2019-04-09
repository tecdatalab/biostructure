import numpy as np
from glumpy import app, gloo, gl, glm
from skimage import measure


import molecule
import reader


class Visualizer():
    
    def __init__(self, myMolecule, labels=None):
        self.vertex = """
                        uniform mat4   model;
                        uniform mat4   view;       // View matrix
                        uniform mat4   projection;   // Projection matrix
                        attribute vec3 position;      // Vertex position
                        attribute vec3 color;
                        varying vec3 v_color;
                        
                        void main()
                        {
                            v_color = color;
                            gl_Position = projection * view * model * vec4(position,1.0);
                        } """
        self.fragment = """
                        varying vec3 v_color;

                        void main()
                        {
                            gl_FragColor = vec4(v_color, 1.0);
                        } """
        self.labels = labels
        self.molecule = myMolecule



    def show(self):

        window = app.Window(512, 512, color=(1,1,1,1))
        phi, theta = 0, 0


        @window.event
        def on_draw(dt):
            nonlocal phi, theta, model
            window.clear()
            gl.glDisable(gl.GL_BLEND)
            gl.glEnable(gl.GL_DEPTH_TEST)
            gl.glEnable(gl.GL_POLYGON_OFFSET_FILL)
            
            points.draw(gl.GL_TRIANGLES, I)

            glm.rotate(model, theta, 1, 0, 0)
            glm.rotate(model, phi, 0, 1, 0)
            
            points['model'] = model


        
        @window.event
        def on_resize(width, height):
            points['projection'] = glm.perspective(45.0, width / float(height), 1.0, 256.0)

        @window.event
        def on_mouse_drag(x, y, dx, dy, buttons):
            nonlocal phi, theta
            print(dx," ",dy)
            phi += dx
            theta += dy
        
            

        Z = self.molecule.data()
        # Marching cubes algorithm to get verts and faces to display
        verts, faces, _,_ = measure.marching_cubes_lewiner(Z)

        verts_index = verts.astype(np.int32)

        if self.labels.any() != None:
            V = np.zeros(len(verts), [("position", np.float32, 3), ("color", np.float32, 3)])
            colors = label2rgb(labels)
            V["color"] = colors[verts_index[:,0], verts_index[:,1], verts_index[:,2]]
        else:
            V = np.zeros(len(verts), [("position", np.float32, 3)])
        
        V["position"] = verts 
        

        V = V.view(gloo.VertexBuffer)
        I = (faces).astype(np.uint32)
        I = I.view(gloo.IndexBuffer)
               
        points = gloo.Program(self.vertex, self.fragment)
        points.bind(V)

        dim = self.molecule.shape()

        model = np.eye(4, dtype=np.float32)
        model = glm.translate(model, -dim[0]/2)
        points['model'] = model

        view = np.eye(4, dtype=np.float32)
        view = glm.translate(view, 0, 0, -dim[0]*2/3)
        points['view'] = view


        app.run()


filename = "../../EMD-1364.map"
myreader = reader.Reader(filename)
myMolecule = myreader.read()

from scipy import ndimage as ndi
from skimage.morphology import watershed, local_maxima
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu
from skimage.color import label2rgb

th = threshold_otsu(myMolecule.data())
binary = myMolecule.data() <= th
distance = ndi.distance_transform_edt(binary)
#local_max = local_maxima(distance , allow_borders=False)
local_max = peak_local_max(distance, indices=False, labels=binary, footprint=np.ones((3,3,3)), exclude_border=1)
markers = ndi.label(local_max, structure=np.ones((3,3,3)))[0] 
labels = watershed(-distance,markers,mask=binary)


v = Visualizer(myMolecule, labels)
v.show()




