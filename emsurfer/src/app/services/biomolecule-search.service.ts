import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Biomolecule } from '../models/biomolecule';
import { BiomoleculeComparison } from '../models/biomolecule-comparison';
import config from '../../config.json';
import { CustomFile } from '../models/custom-file';

@Injectable({
  providedIn: 'root'
})
export class BiomoleculeSearchService {
  constructor(private httpClient: HttpClient) {}

  readonly API_URL = config.api_url;

  getBiomolecule(emdbId: number): Promise<void | Biomolecule> {
    return this.httpClient
      .get(this.API_URL + '/search/' + emdbId)
      .toPromise()
      .then(response => {
        const object = response as Biomolecule;
        object.image_url = this.API_URL + object.image_url;
        return object;
      })
      .catch(this.handleError);
  }

  getZernikeDescriptors(
    emdbId: number,
    contourRepresentation: number
  ): Promise<any> {
    return this.httpClient
      .get(
        this.API_URL + '/search/zernike/' + emdbId + '/' + contourRepresentation
      )
      .toPromise()
      .then()
      .catch(this.handleError);
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
  ): Promise<any> {
    return this.httpClient
      .get(
        this.API_URL +
          '/search/' +
          emdbId +
          '/' +
          isVolumeFilterOn +
          '/' +
          minRes +
          '/' +
          maxRes
      )
      .toPromise()
      .then((data: any) => {
        for (let item of data.results) {
          item.biomolecule.image_url =
            this.API_URL + item.biomolecule.image_url;
        }
        console.log(data);
        return data;
      })
      .catch(this.handleError);
  }

  private handleError(error: any) {
    let errMsg = error.message
      ? error.message
      : error.status
      ? `${error.status} - ${error.statusText}`
      : 'Server error';
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
