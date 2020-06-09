from process_mrc.miscellaneou import get_mrc_segments
from process_graph.process_segment_faces import get_n_points_cube
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TKAgg')


segments = get_mrc_segments("../../maps/1010/EMD-1010.map", 7, 3, 1)
face_points = get_n_points_cube(segments[0].mask, 3, 0)

x, y, z = [],[],[]
for i in face_points:
    x.append(i[0])
    y.append(i[1])
    z.append(i[2])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c='b', marker='o')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()