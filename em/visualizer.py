import numpy as np
from glumpy import app, gloo, gl, glm
from glumpy.transforms import Trackball, Position
from glumpy.ext import png
from skimage import measure
from skimage.color import label2rgb
from biopandas.pdb import PandasPdb 
from skimage.filters import threshold_otsu
from skimage.measure import regionprops
from collections import Counter



import reader
import processing

class Visualizer():
    
    def __init__(self, myMolecule, level=None):
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
                        gl_FragColor = vec4(color_gamma, 0.8);
                    }
                    """
        self.molecule = myMolecule
        data = self.molecule.data()

        self.th_level = np.min(data) + level * (np.max(data) - np.min(data)) if level is not None else threshold_otsu(data)
        verts, faces, _,_ = measure.marching_cubes_lewiner(data, self.th_level)

        self.verts = verts
        self.faces = faces      
        self.atoms = None
        self.labels = None


    def show(self, export=False, start_angle=None, time=1):

        faces = self.faces
        verts = self.verts
        verts_index = verts.astype(np.int32)

        if self.labels.any() != None:
            V = np.zeros(len(verts), [("position", np.float32, 3),  ("normal", np.float32, 3), ("color", np.float32, 3)])
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

        I = (faces).astype(np.uint32)

        ###ATOMS CODE
        if self.atoms is not None:
            atoms_coords = self.atoms
            amin, amax =  np.min(atoms_coords["position"]), np.max(atoms_coords["position"])
            atoms_norm = 2*(atoms_coords["position"]-amin)/(amax-amin) - 1
            ####Generalte atoms vertices
            
            slices = 32 + 1
            stacks = 32 + 1
            n = slices*stacks
            V_atoms = np.zeros(len(atoms_norm)*n, [("position", np.float32, 3),  ("normal", np.float32, 3), ("color", np.float32, 3)])
            sphere_verts = np.zeros(n, [("position", np.float32, 3), ("normal", np.float32, 3), ("color", np.float32, 3)])
            radius = 0.02
            theta1 = np.repeat(np.linspace(0,     np.pi, stacks, endpoint=True), slices)
            theta2 = np.tile(np.linspace(0, 2 * np.pi, slices, endpoint=True), stacks)
            sphere_verts["position"][:,0] =  np.sin(theta1) * np.cos(theta2) * radius
            sphere_verts["position"][:,1] =  np.cos(theta1) * radius
            sphere_verts["position"][:,2] =  np.sin(theta1) * np.sin(theta2) * radius
            sphere_verts = np.repeat(sphere_verts[None,:], len(atoms_norm), axis=0).reshape(-1)
            atoms_norm = np.repeat(atoms_norm, n, axis=0)
            V_atoms["position"][:,0] = atoms_norm[:,2] + sphere_verts["position"][:,2]
            V_atoms["position"][:,1] = atoms_norm[:,1] + sphere_verts["position"][:,1]
            V_atoms["position"][:,2] = atoms_norm[:,0] + sphere_verts["position"][:,0]
            V_atoms["normal"] = V_atoms["position"]
            V_atoms["color"] = np.repeat(np.array([[0.0, 0.1, 1.0]]),len(atoms_norm),axis=0)

            vamin, vamax =  V_atoms["position"].min(), V_atoms["position"].max()
            V_atoms["position"] = 2*(V_atoms["position"]-vamin)/(vamax-vamin) - 1
            
            indices = []
            initial_offset = 0
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

            V = np.concatenate((V_atoms, V), axis=0)
            I = np.concatenate((indices,(faces).astype(np.uint32)+len(V_atoms)), axis=0)
        
        
        V = V.view(gloo.VertexBuffer)
        I = I.view(gloo.IndexBuffer)

        points = gloo.Program(self.vertex, self.fragment)
        points.bind(V)
        trackball = Trackball(Position("position"))
        points['transform'] = trackball
        trackball.zoom = 15

        if start_angle != None:
            trackball.theta, trackball.phi = start_angle[0], start_angle[1]
        else:
            trackball.theta, trackball.phi = 0, 0

        
        window = app.Window(width=512, height=512, color=(1,1,1,1))

        framebuffer = np.zeros((window.height, window.width * 3), dtype=np.uint8)

        def update():
            model = points['transform']['model'].reshape(4,4)
            view  = points['transform']['view'].reshape(4,4)
            points['m_view']  = view
            points['m_model'] = model
            points['m_normal'] = np.array(np.matrix(np.dot(view, model)).I.T)

        @window.event
        def on_mouse_drag(x, y, dx, dy, button):
            update()
            
        @window.event
        def on_draw(dt):
            window.clear()
            points.draw(gl.GL_TRIANGLES, I)
            
            
        @window.event
        def on_init():
            gl.glEnable(gl.GL_DEPTH_TEST)
            update()

        if export:
            frame_count = 1

            @window.timer(0.1) 
            def timer(elapsed):
                nonlocal frame_count
                points['transform'].phi +=2
                update()
                gl.glReadPixels(0, 0, window.width, window.height, gl.GL_RGB, gl.GL_UNSIGNED_BYTE, framebuffer)
                png.from_array(framebuffer, 'RGB').save('export/frame'+str(frame_count)+'.png')
                frame_count+=1

                
        window.attach(points['transform'])
        if export:
            app.run(duration=time)
        else:
            app.run()


    def add_structure(self, filename):
        ppdb = PandasPdb()
        ppdb.read_pdb(filename)
        atoms = ppdb.df["ATOM"][ppdb.df['ATOM']['atom_name'] == 'CA']
        self.atoms = np.zeros(len(atoms),  [("id", np.int32, 1), ("position", np.float32, 3)])
        self.atoms["id"] = atoms["atom_number"].values
        self.atoms["position"] = atoms[["z_coord", "y_coord", "x_coord"]].values

    def show_atom_correlation(self):
        if self.atoms is None:
            raise ValueError("Atomic structure is not present")
        if self.labels is None:
            raise ValueError("Need to segmentate map first")
        else:
            '''
            molecule = self.molecule
            atoms_coords = self.atoms["position"]
            verts = self.verts
            zlen,ylen,xlen = molecule.cell_dim()
            mz, my, mx = molecule.grid_size()
            map_size = np.array([mz, my, mx])
            cell_dim = np.array([xlen,ylen,zlen])
            voxel_len = cell_dim/map_size
            map_min = np.min(verts, axis=0)
            map_max = np.max(verts, axis=0)
            voxel_padding = map_size - (map_max-map_min)
            padding_adjustment = voxel_padding * voxel_len / 2
            pdb_min = np.min(atoms_coords, axis=0)
            pdb_max = np.max(atoms_coords, axis=0)
            adjustment = (pdb_min - padding_adjustment)/voxel_len
            '''
            verts = self.verts
            atoms_coords = self.atoms["position"]
            map_min = np.min(verts, axis=0)
            map_max = np.max(verts, axis=0)
            pdb_min = np.min(atoms_coords, axis=0)
            pdb_max = np.max(atoms_coords, axis=0)
            
            atoms_labels = dict()

            for atom in self.atoms:
                atoms_labels[atom["id"]] = []
                for region in regionprops(self.labels):
                    min_z,min_y,min_x,max_z,max_y,max_x= region.bbox
                    #min_box = np.array([min_z,min_y,min_x]) + adjustment
                    #max_box = np.array([max_z,max_y,max_x]) + adjustment
                    #coord = atom["position"]
                    min_box = 2*(np.array([min_z,min_y,min_x])-map_min)/(map_max-map_min) - 1
                    max_box = 2*(np.array([max_z,max_y,max_x])-map_min)/(map_max-map_min) - 1
                    coord = 2*(atom["position"]-pdb_min)/(pdb_max-pdb_min) -1
                    if np.all(coord <= max_box) and np.all(coord >= min_box):
                        atoms_labels[atom["id"]].append(region.label)
            atoms_labels = self.assign_atom_knn(atoms_labels, 50)
            print(Counter(atoms_labels.values()))
           
    def assign_atom_knn(self, atoms_dict, k):
        atomic_structure = self.atoms
        new_atoms_dict = dict()
        atoms_distance = np.zeros(len(atomic_structure), [("id", np.int32, 1), ("distance", np.float32, 1)])
        atoms_distance["id"] = atomic_structure["id"][:]
        for key in atoms_dict:
            distance = lambda x: np.linalg.norm(atomic_structure[atomic_structure["id"]==key]['position']-x, axis=1)
            atoms_distance["distance"] = distance(atomic_structure["position"])
            sorted_distance = np.sort(atoms_distance, order=["distance"])
            selected_keys = sorted_distance["id"][1:k+1]
            selected_labels = np.array([label for label_list in [atoms_dict[selected_key] for selected_key in selected_keys] for label in label_list])
            if selected_labels != []:
                label = np.bincount(selected_labels).argmax()
            else:
                label = 0
            new_atoms_dict[key] = label
        return new_atoms_dict
            


    def segmentate(self, step_sigma, steps):
        molecule = self.molecule
        step_maps = processing.gaussian_step_filtering(molecule, step_sigma, steps)
        self.labels = processing.watershed_segmentation(molecule, self.th_level, step_maps, steps)








#Read molecule map from file
mapReader = reader.Reader()
#Open file
mapReader.open("../maps/1364/EMD-1364.map")
#Get map object
myMap = mapReader.read()
# Create visualizer with a map surface threshold level
# Otherwise use otsu threshold
v= Visualizer(myMap)
#v = Visualizer(myMap)
# Watershed 
v.segmentate(step_sigma=0.3, steps=3)
# add corresponding atomic structure
v.add_structure("../maps/1364/pdb1pn6.ent")
v.show()
v.show_atom_correlation()


