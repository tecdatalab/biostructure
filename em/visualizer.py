import numpy as np
from glumpy import app, gloo, gl, glm
from skimage import measure
from skimage.color import label2rgb


import molecule
import reader


class Visualizer():
    
    def __init__(self, myMolecule, labels=None):
        self.vertex = """
                        uniform mat4      u_model;
                        uniform mat4      u_view;       // View matrix
                        uniform mat4      u_projection;   // Projection matrix

                        attribute vec3      a_position;      // Vertex position
                        attribute vec3      a_normal;
                        attribute vec3      u_color;

                        varying vec3      v_color;
                        varying vec3      v_normal;
                        varying vec3      v_position;

                        
                        
                        void main()
                        {
                            v_color = u_color;
                            v_normal   = a_normal;
                            v_position = a_position;

                            gl_Position = u_projection * u_view * u_model * vec4(a_position,1.0);
                        } """

        self.fragment = """
                        uniform mat4      u_model;           // Model matrix
                        uniform mat4      u_view;            // View matrix
                        uniform mat4      u_normal;          // Normal matrix
                        uniform mat4      u_projection;      // Projection matrix
                        uniform vec3      u_light_position;  // Light position
                        uniform vec3      u_light_intensity; // Light intensity
                        
                        varying vec3      v_normal;          // Interpolated normal (in)
                        varying vec3      v_position;        // Interpolated position (in)
                        varying vec3      v_color;           // Interpolated color (in)

                        void main()
                        {
                            // Calculate normal in world coordinates
                            vec3 normal = normalize(u_normal * vec4(v_normal,1.0)).xyz;
                            // Calculate the location of this fragment (pixel) in world coordinates
                            vec3 position = vec3(u_view*u_model * vec4(v_position, 1.0));
                            // Calculate the vector from this pixels surface to the light source
                            vec3 surfaceToLight = u_light_position - position;
                            // Calculate the cosine of the angle of incidence (brightness)
                            float brightness = dot(normal, surfaceToLight) / (length(surfaceToLight) * length(normal));
                            brightness = max(min(brightness,1.0),0.0);

                            // Final color
                            //gl_FragColor = vec4(v_color, 1.0) * (0.1 + 0.9*brightness * vec4(u_light_intensity, 1));
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

            phi = 0
            theta = 0

            view = points['u_view'].reshape(4,4)            
            points['u_model'] = model
            points['u_normal'] = np.array(np.matrix(np.dot(view, model)).I.T)


        
        @window.event
        def on_resize(width, height):
            points['u_projection'] = glm.perspective(45.0, width / float(height), 1.0, 256.0)

        @window.event
        def on_mouse_drag(x, y, dx, dy, buttons):
            nonlocal phi, theta
            phi += dx
            theta += dy
        
            

        Z = self.molecule.data()
        # Marching cubes algorithm to get verts and faces to display
        verts, faces, _,_ = measure.marching_cubes_lewiner(Z)

        verts_index = verts.astype(np.int32)

        if self.labels.any() != None:
            V = np.zeros(len(verts), [("a_position", np.float32, 3), ("u_color", np.float32, 3)])
            colors = label2rgb(self.labels)
            V["u_color"] = colors[verts_index[:,0], verts_index[:,1], verts_index[:,2]]
        else:
            V = np.zeros(len(verts), [("a_position", np.float32, 3)])
        
        V["a_position"] = verts 
        

        V = V.view(gloo.VertexBuffer)
        I = (faces).astype(np.uint32)
        I = I.view(gloo.IndexBuffer)
               
        points = gloo.Program(self.vertex, self.fragment)
        points.bind(V)

        dim = self.molecule.shape()

        model = np.eye(4, dtype=np.float32)
        model = glm.translate(model, -dim[0]/2)
        points['u_model'] = model

        view = np.eye(4, dtype=np.float32)
        view = glm.translate(view, 0, 0, -dim[0]*4/5)
        points['u_view'] = view

        points["u_light_position"] = -2,-2,2
        points["u_light_intensity"] = 1,1,1


        app.run()





