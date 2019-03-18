import { Injectable } from '@angular/core';
import { Biomolecule } from '../models/biomolecule';
import { BiomoleculeComparison } from '../models/biomolecule-comparison';

@Injectable({
  providedIn: 'root'
})
export class BiomoleculeSearchService {
  constructor() {}

  getBiomolecule(emdbId: number) {
    const newBiomolecule = new Biomolecule();
    newBiomolecule.emdb_id = emdbId;
    newBiomolecule.img_url = '../../../assets/img/test_img.gif';
    newBiomolecule.name = 'Lorem ipsum et doloren';
    return newBiomolecule;
  }

  getZernikeDescriptors(emdbId: number, contourRepresentation: number) {
    return [1, 2, 3, 4, 3, 10, 0];
  }

  getSimilarBioMolecules(
    emdbId: number,
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
      newBiomolecule.biomolecule.name = 'Lorem ipsum et doloren';
      newBiomolecule.euc_distance = 1;
      newBiomolecule.ratio_of_volume = 0.5;
      newBiomolecule.resolution = 10.1;
      results.push(newBiomolecule);
    }
    return results;
  }
}
