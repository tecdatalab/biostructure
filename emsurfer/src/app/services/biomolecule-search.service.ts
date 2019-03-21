import { Injectable } from '@angular/core';
import { Biomolecule } from '../models/biomolecule';
import { BiomoleculeComparison } from '../models/biomolecule-comparison';
import { CustomFile } from '../models/custom-file';

@Injectable({
  providedIn: 'root'
})
export class BiomoleculeSearchService {
  constructor() {}

  getBiomolecule(emdbId: number) {
    const newBiomolecule = new Biomolecule();
    newBiomolecule.emdb_id = emdbId;
    newBiomolecule.img_url = '../../../assets/img/test_img.gif';
    newBiomolecule.pdb_url = 'http://www.ebi.ac.uk/pdbe/entry/emdb/EMD-1413';
    newBiomolecule.name = 'Lorem ipsum et doloren';
    return newBiomolecule;
  }

  getZernikeDescriptors(emdbId: number, contourRepresentationId: number) {
    if (emdbId === 5555) {
      return [10, 11, 0, 3, 2, 4, 5, 15];
    }
    return [1, 2, 3, 4, 3, 10, 0];
  }

  getZernikeDescriptorsByMapId(
    emMapId: number,
    contourLevel: number,
    contourRepresentationId: number
  ) {
    return [1, 2, 1, 2, 1, 2, 3];
  }

  getSimilarBioMolecules(
    emdbId: number,
    contourRepresentationId: number,
    isVolumeFilterOn: boolean,
    minRes: number,
    maxRes: number
  ) {
    const results = [];
    for (let i = 0; i < 5; i++) {
      const newBiomolecule = new BiomoleculeComparison();
      newBiomolecule.biomolecule = new Biomolecule();
      newBiomolecule.biomolecule.emdb_id = emdbId;
      newBiomolecule.biomolecule.img_url = '../../../assets/img/test_img.gif';
      newBiomolecule.biomolecule.pdb_url =
        'http://www.ebi.ac.uk/pdbe/entry/emdb/EMD-1413';
      newBiomolecule.biomolecule.name = 'Lorem ipsum et doloren';
      newBiomolecule.euc_distance = 1;
      newBiomolecule.ratio_of_volume = 0.5;
      newBiomolecule.resolution = 10.1;
      results.push(newBiomolecule);
    }
    return results;
  }

  getBatchBiomolecules(
    fileId: number,
    contourRepresentationId: number,
    isVolumeFilterOn: boolean,
    topResults: number
  ) {
    const files = [];
    for (let i = 0; i < 5; i++) {
      const f = new CustomFile();
      f.filename = 'EMDB-' + i + i + i + i + '.hit';
      f.path = 'assets/test_files/test_result.hit';
      files.push(f);
    }
    return files;
  }

  getBatchBiomoleculesByFileId(
    fileId: number,
    contourRepresentationId: number,
    isVolumeFilterOn: boolean,
    topResults: number
  ) {
    const files = [];
    for (let i = 0; i < 5; i++) {
      const f = new CustomFile();
      f.filename = 'EMDBF-' + i + i + i + i + '.hit';
      f.path = 'assets/test_files/test_result.hit';
      files.push(f);
    }
    return files;
  }
}
