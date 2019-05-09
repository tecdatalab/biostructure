import numpy as np
from glumpy import app, gloo, gl, glm
from skimage import measure
from skimage.color import label2rgb
from biopandas.pdb import PandasPdb 


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
                        attribute vec3      a_color;

                        varying vec3      v_color;
                        varying vec3      v_normal;
                        varying vec3      v_position;

                        
                        
                        void main()
                        {
                            v_color = a_color;
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
                            brightness = max(min(brightness,1.0),0.25);

                            // Final color
                            gl_FragColor = vec4(v_color, 1.0) * (0.1 + 0.9*brightness * vec4(u_light_intensity, 1));
                            //gl_FragColor = vec4(v_color, 1.0);
                        } """

        self.labels = labels
        self.molecule = myMolecule
        self.atoms = None



    def show(self):

        window = app.Window(512, 512, color=(1,1,1,1))
        phi, theta = 0, 0


        @window.event
        def on_draw(dt):
            nonlocal phi, theta, model
            window.clear()
            gl.glDisable(gl.GL_BLEND)
            gl.glEnable(gl.GL_DEPTH_TEST)
            #gl.glEnable(gl.GL_POLYGON_OFFSET_FILL)
            
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
            V = np.zeros(len(verts), [("a_position", np.float32, 3),  ("a_normal", np.float32, 3), ("a_color", np.float32, 3)])
            colors = label2rgb(self.labels[verts_index[:,0], verts_index[:,1], verts_index[:,2]])
            V["a_color"] = colors
        else:
            V = np.zeros(len(verts), [("a_position", np.float32, 3),  ("a_normal", np.float32, 3)])
        
        V["a_position"] = verts 

        triangles = verts[faces]
        n = np.cross( triangles[::,1 ] - triangles[::,0] , triangles[::,2 ] - triangles[::,0] )
        n = n/np.linalg.norm(n, axis=1)[:,None]
        V["a_normal"][ faces[:,0] ] += n 
        V["a_normal"][ faces[:,1] ] += n 
        V["a_normal"][ faces[:,2] ] += n 
        V["a_normal"] = V["a_normal"] / np.linalg.norm(V["a_normal"], axis=1)[:,None]
   

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
        view = glm.translate(view, 0, 0, -dim[2]*4/5)
        points['u_view'] = view

        points["u_light_position"] = 0,-1000,0
        points["u_light_intensity"] = 1,1,1


        app.run()
        return verts

    def add_structure(self, filename):
        ppdb =PandasPdb()
        ppdb.read_pdb(filename)
        atoms = ppdb.df["ATOM"][ppdb.df['ATOM']['atom_name'] == 'CA']
        self.atoms = atoms[["z_coord", "y_coord", "x_coord"]].values
        return self.atoms



r = reader.Reader("../../EMD-1364.map")
m = r.read()
#cube = np.zeros((130,130,130))
#cube[50:80,50:80,50:80] = 1
#m.set_data(cube)
import processing
from skimage.measure import regionprops
s = processing.gaussian_smooth(m,1)
l = processing.watershed_segmentation(m)
v = Visualizer(m,l)
atoms_coords = v.add_structure("../../pdb1mi6.ent")
zlen,ylen,xlen = m.cell_dim()
mz, my, mx = m.grid_size()
map_size = np.array([mz, my, mx])
cell_dim = np.array([xlen,ylen,zlen])
voxel_len = cell_dim/map_size
print(voxel_len)
print(cell_dim)


verts =v.show()
map_min = np.min(verts, axis=0)
map_max = np.max(verts, axis=0)
voxel_padding = map_size - (map_max-map_min)
padding_adjustment = voxel_padding * voxel_len / 2
print(padding_adjustment)
pdb_min = np.min(atoms_coords, axis=0)
adjustment = (pdb_min - padding_adjustment)/voxel_len
print(adjustment)

atoms_found=1
for coord in atoms_coords:
    for region in regionprops(l):
        min_z,min_y,min_x,max_z,max_y,max_x= region.bbox
        min_box = np.array([min_z,min_y,min_x]) + adjustment
        max_box = np.array([max_z,max_y,max_x]) + adjustment
        if np.all(coord <= max_box) and np.all(coord >= min_box):
            print("Found atom %d with coords (%f, %f, %f) in region %d " % (atoms_found, coord[0], coord[1], coord[2], region.label))
            atoms_found+=1


