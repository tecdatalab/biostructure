import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Biomolecule } from '../models/biomolecule';
import { BiomoleculeComparison } from '../models/biomolecule-comparison';
import config from '../../config.json';

@Injectable({
  providedIn: 'root'
})
export class BiomoleculeSearchService {
  constructor(private httpClient: HttpClient) {}

  readonly API_URL = config.api_url;

  getBiomolecule(emdbId: number): Promise<void | Biomolecule> {
    /*
    const newBiomolecule = new Biomolecule();
    newBiomolecule.emdb_id = emdbId;
    newBiomolecule.img_url = '../../../assets/img/test_img.gif';
    newBiomolecule.name = 'Lorem ipsum et doloren';
    */
    return this.httpClient
      .get(this.API_URL + '/search/' + emdbId)
      .toPromise()
      .then(response => {
        let object = response as Biomolecule;
        object.image_url = this.API_URL + object.image_url;
        return object;
      })
      .catch(this.handleError);
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
      newBiomolecule.biomolecule.id = emdbId;
      newBiomolecule.biomolecule.image_url = '../../../assets/img/test_img.gif';
      newBiomolecule.biomolecule.acronym = 'Lorem ipsum et doloren';
      newBiomolecule.euc_distance = 1;
      newBiomolecule.ratio_of_volume = 0.5;
      newBiomolecule.resolution = 10.1;
      results.push(newBiomolecule);
    }
    return results;
  }

  private handleError(error: any) {
    let errMsg = error.message
      ? error.message
      : error.status
      ? `${error.status} - ${error.statusText}`
      : 'Server error';
    console.error('ERRRRRRRRRRROR ' + errMsg);
  }
}
