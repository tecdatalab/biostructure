import numpy as np
import mrcfile
import os
from shutil import copyfile

from general_utils.download_utils import download_pdb
from general_utils.pdb_utils import get_chains_pdb
from to_mrc.pdb_2_mrc import pdb_to_mrc_chains


def get_point_draw(shape, pos, point_percentage):
  long_x = int((shape[0] * point_percentage) / 2)
  long_y = int((shape[1] * point_percentage) / 2)
  long_z = int((shape[2] * point_percentage) / 2)

  print(long_x, long_y, long_z)

  x_range = [max(pos[0] - long_x, 0), min(shape[0] - 1, pos[0] + long_x)]
  y_range = [max(pos[1] - long_y, 0), min(shape[1] - 1, pos[1] + long_y)]
  z_range = [max(pos[2] - long_z, 0), min(shape[2] - 1, pos[2] + long_z)]

  return (x_range, y_range, z_range)


def create_ani_expe_1a(dir_path, name_pdb, chains_father,
                       chains_child, match_list, point_percentage, center_father,
                       center_child, test_center_father, test_center_child, resolution,
                       is_mrc_files_existed=False):
  if not is_mrc_files_existed:
    if not os.path.isdir(dir_path):
      os.mkdir(dir_path)

    download_pdb('{1}', '{0}/{1}.pdb'.format(dir_path, name_pdb))

    # Maps creation
    chains = get_chains_pdb('{0}/{1}.pdb'.format(dir_path, name_pdb))
    pdb_to_mrc_chains(True, False, resolution, '{0}/{1}.pdb'.format(dir_path, name_pdb), dir_path, chains, len(chains))
    os.remove('{0}/{1}.pdb'.format(dir_path, name_pdb))

  work_path = os.path.abspath(dir_path) + '/' + name_pdb

  original_mrc = mrcfile.open('{}/{}.mrc'.format(work_path, name_pdb))
  shape_original = original_mrc.data.shape

  create_point(center_father, name_pdb, original_mrc, point_percentage, shape_original, work_path, 'father_pc')
  create_point(test_center_father, name_pdb, original_mrc, point_percentage, shape_original, work_path, 'father_pt')
  create_point(center_child, name_pdb, original_mrc, point_percentage, shape_original, work_path, 'child_pc')
  create_point(test_center_child, name_pdb, original_mrc, point_percentage, shape_original, work_path, 'child_pt')

  #create pymol file father
  if os.path.exists('{}/script_father.pml'.format(work_path)):
    os.remove('{}/script_father.pml'.format(work_path))

  f = open('{}/script_father.pml'.format(work_path), 'a+')

  #Load mrc files
  f.write('load {}/{}.mrc'.format(work_path, name_pdb))
  f.write('\n')
  for i in chains_father:
    f.write('load {}/{}_{}.mrc'.format(work_path, name_pdb, i))
    f.write('\n')
  f.write('load {}/{}_father_pc.mrc'.format(work_path, name_pdb))
  f.write('\n')
  f.write('load {}/{}_father_pt.mrc'.format(work_path, name_pdb))
  f.write('\n\n\n')

  #Gen volume
  f.write('volume {0}_volume, {0}\n'.format(name_pdb))
  f.write('show volume, {0}_volume\n'.format(name_pdb))
  f.write('\n\n\n')

  #Gen surface
  for i in chains_father:
    f.write('isosurface {0}_{1}_surface, {0}_{1}'.format(name_pdb, i))
    f.write('\n')
    f.write('hide surface, {0}_{1}_surface'.format(name_pdb, i))
    f.write('\n')

  f.write('isosurface {0}_father_pc_surface, {0}_father_pc'.format(name_pdb))
  f.write('\n')
  f.write('color red, {0}_father_pc_surface'.format(name_pdb))
  f.write('\n')
  f.write('isosurface {0}_father_pt_surface, {0}_father_pt'.format(name_pdb))
  f.write('\n')
  f.write('color white, {0}_father_pt_surface'.format(name_pdb))
  f.write('\n\n\n')

  #Gen Images
  f.write('png {0}/{1}.png'.format(work_path, 1))
  f.write('\n\n\n')



def create_point(center_father, name_pdb, original_mrc, point_percentage, shape_original, work_path, file_name):
  copyfile('{}/{}.mrc'.format(work_path, name_pdb), '{}/{}_{}.mrc'.format(work_path, name_pdb, file_name))
  with mrcfile.open('{}/{}_{}.mrc'.format(work_path, name_pdb, file_name), mode='r+') as mrc:
    # mrc.set_data(np.zeros(shape_original, dtype=np.int8))
    mrc.data[:, :, :] = 0
    x_range, y_range, z_range = get_point_draw(shape_original, center_father, point_percentage)
    mrc.data[x_range[0]:x_range[1], y_range[0]:y_range[1], z_range[0]:z_range[1]] = np.max(original_mrc.data)
