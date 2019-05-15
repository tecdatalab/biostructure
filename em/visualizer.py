import numpy as np
from glumpy import app, gloo, gl, glm
from glumpy.transforms import Trackball, Position
from skimage import measure
from skimage.color import label2rgb
from biopandas.pdb import PandasPdb 


import reader

class Visualizer():
    
    def __init__(self, myMolecule, labels=None):
        self.vertex = """
                    uniform mat4 m_model;
                    uniform mat4 m_view;
                    uniform mat4 m_normal;
                    attribute vec3 position;
                    attribute vec3 normal;
                    attribute vec3 color;
                    varying vec3 v_normal;
                    varying vec3 v_position;
                    varying vec3 v_color;
                    void main()
                    {

                        gl_Position = <transform>;
                        vec4 P = m_view * m_model* vec4(position, 1.0);
                        v_position = P.xyz / P.w;
                        v_normal = vec3(m_normal * vec4(normal,0.0));
                        v_color = color;
                    }
                    """

        self.fragment = """
                    varying vec3 v_normal;
                    varying vec3 v_position;
                    varying vec3 v_color;
                    const vec3 light_position = vec3(1.0,1.0,1.0);
                    vec3 ambient_color = v_color;
                    const vec3 diffuse_color = vec3(0.125, 0.125, 0.125);
                    const vec3 specular_color = vec3(1.0, 1.0, 1.0);
                    const float shininess = 128.0;
                    const float gamma = 2.2;
                    void main()
                    {
                        vec3 normal= normalize(v_normal);
                        vec3 light_direction = normalize(light_position + v_position);
                        float lambertian = max(dot(light_direction,normal), 0.0);
                        float specular = 0.0;
                        if (lambertian > 0.0)
                        {
                            vec3 view_direction = normalize(v_position);
                            vec3 half_direction = normalize(light_direction + view_direction);
                            float specular_angle = max(dot(half_direction, normal), 0.0);
                            specular = pow(specular_angle, shininess);
                        }
                        vec3 color_linear = ambient_color +
                                            lambertian * diffuse_color +
                                            specular * specular_color;
                        vec3 color_gamma = pow(color_linear, vec3(1.0/gamma));
                        gl_FragColor = vec4(color_gamma, 0.7);
                    }
                    """
        self.labels = labels
        self.molecule = myMolecule
        self.atoms = None


    def show(self):
        data = self.molecule.data()
        verts, faces, _,_ = measure.marching_cubes_lewiner(data)

        verts_index = verts.astype(np.int32)

        if self.labels.any() != None:
            V = np.zeros(len(verts), [("position", np.float32, 3),  ("normal", np.float32, 3), ("color", np.float32, 3)])
            print(self.labels)
            colors = label2rgb(self.labels[verts_index[:,0], verts_index[:,1], verts_index[:,2]])
            V["color"] = colors * 0.5
        else:
            V = np.zeros(len(verts), [("position", np.float32, 3),  ("normal", np.float32, 3)])
        
        T = verts[faces]
        N = np.cross(T[::,1 ]-T[::,0], T[::,2]-T[::,0])
        L = np.sqrt(N[:,0]**2+N[:,1]**2+N[:,2]**2)
        N /= L[:, np.newaxis]
        normals = np.zeros(verts.shape)
        normals[faces[:,0]] += N
        normals[faces[:,1]] += N
        normals[faces[:,2]] += N
        L = np.sqrt(normals[:,0]**2+normals[:,1]**2+normals[:,2]**2)
        normals /= L[:, np.newaxis]

        vmin, vmax =  verts.min(), verts.max()
        verts_norm = 2*(verts-vmin)/(vmax-vmin) - 1

        V["position"] = verts_norm
        V["normal"] = normals / np.linalg.norm(normals, axis=1)[:,None]

        ###ATOMS CODE
        atoms_coords = self.add_structure("../../pdb1mi6.ent")
        amin, amax =  atoms_coords.min(), atoms_coords.max()
        atoms_norm = 2*(atoms_coords-amin)/(amax-amin) - 1
        ####Generalte atoms vertices
        
        slices = 8 + 1
        stacks = 8 + 1
        n = slices*stacks
        V_atoms = np.zeros(len(atoms_norm)*n, [("position", np.float32, 3),  ("normal", np.float32, 3), ("color", np.float32, 3)])
        sphere_verts = np.zeros(n, [("position", np.float32, 3), ("normal", np.float32, 3), ("color", np.float32, 3)])
        radius = 0.03
        theta1 = np.repeat(np.linspace(0,     np.pi, stacks, endpoint=True), slices)
        theta2 = np.tile(np.linspace(0, 2 * np.pi, slices, endpoint=True), stacks)
        sphere_verts["position"][:,0] =  np.sin(theta1) * np.cos(theta2) * radius
        sphere_verts["position"][:,1] =  np.cos(theta1) * radius
        sphere_verts["position"][:,2] =  np.sin(theta1) * np.sin(theta2) * radius
        print (sphere_verts.shape)
        sphere_verts = np.repeat(sphere_verts[None,:], len(atoms_norm), axis=0).reshape(-1)
        print (sphere_verts[0])
        print (sphere_verts.shape)
        print (atoms_norm.shape)
        atoms_norm = np.repeat(atoms_norm, n, axis=0)
        print (atoms_norm.shape)
        V_atoms["position"][:,0] = atoms_norm[:,2] + sphere_verts["position"][:,2]
        V_atoms["position"][:,1] = atoms_norm[:,1] + sphere_verts["position"][:,1]
        V_atoms["position"][:,2] = atoms_norm[:,0] + sphere_verts["position"][:,0]
        V_atoms["normal"] = V_atoms["position"]
        V_atoms["color"] = np.repeat(np.array([[0.0, 0.1, 1.0]]),len(atoms_norm),axis=0)

        vamin, vamax =  V_atoms["position"].min(), V_atoms["position"].max()
        V_atoms["position"] = 2*(V_atoms["position"]-vamin)/(vamax-vamin) - 1
        
        V = np.concatenate((V, V_atoms), axis=0)

        indices = []
        initial_offset = len(verts)
        for atom in range(len(atoms_coords)):
            for i in range(stacks-1):
                for j in range(slices-1):
                    indices.append(i*(slices) + j + initial_offset)
                    indices.append(i*(slices) + j+1 + initial_offset)
                    indices.append(i*(slices) + j+slices+1 + initial_offset)
                    indices.append(i*(slices) + j+slices + initial_offset)
                    indices.append(i*(slices) + j+slices+1+ initial_offset)
                    indices.append(i*(slices) + j+ initial_offset)
            initial_offset+=n
        indices = np.array(indices, dtype=np.uint32).reshape(-1,3)

        I = np.concatenate(((faces).astype(np.uint32),indices), axis=0)
        
        
        V = V.view(gloo.VertexBuffer)
        #I = (faces).astype(np.uint32)
        I = I.view(gloo.IndexBuffer)

        points = gloo.Program(self.vertex, self.fragment)
        points.bind(V)
        trackball = Trackball(Position("position"))
        points['transform'] = trackball
        trackball.theta, trackball.phi, trackball.zoom = 0, 0, 15

        window = app.Window(width=1024, height=1024, color=(1,1,1,1))

        def update():
            model = points['transform']['model'].reshape(4,4)
            view  = points['transform']['view'].reshape(4,4)
            points['m_view']  = view
            points['m_model'] = model
            points['m_normal'] = np.array(np.matrix(np.dot(view, model)).I.T)
            
        @window.event
        def on_draw(dt):
            window.clear()
            points.draw(gl.GL_TRIANGLES, I)

        @window.event
        def on_mouse_drag(x, y, dx, dy, button):
            update()
            
        @window.event
        def on_init():
            gl.glEnable(gl.GL_DEPTH_TEST)
            update()

        window.attach(points['transform'])
        app.run()
        return verts

    def add_structure(self, filename):
        ppdb =PandasPdb()
        ppdb.read_pdb(filename)
        atoms = ppdb.df["ATOM"][ppdb.df['ATOM']['atom_name'] == 'CA']
        self.atoms = atoms[["z_coord", "y_coord", "x_coord"]].values
        return self.atoms



r = reader.Reader("../../EMD-1010.map")
m = r.read()
#cube = np.zeros((130,130,130))
#cube[50:80,50:80,50:80] = 1
#m.set_data(cube)
import processing
from skimage.measure import regionprops
s = processing.gaussian_smooth(m,1.5)
l = processing.watershed_segmentation(s)
v = Visualizer(m,l)
v.show()
'''
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
'''

