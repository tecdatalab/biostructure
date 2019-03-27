import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import { Biomolecule } from "../models/biomolecule";
import { BenchmarkResult } from "../models/benchmark-result";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class BiomoleculeSearchService {
  constructor(private httpClient: HttpClient) {}

  readonly API_URL = config.api_url;

  getBiomolecule(emdbId: number): Promise<void | Biomolecule> {
    return this.httpClient
      .get(this.API_URL + "/search/" + emdbId)
      .toPromise()
      .then((response: Biomolecule) => {
        //const object = response as Biomolecule;
        console.log(response);
        response.image_url = this.API_URL + response.image_url;
        return response;
      })
      .catch(this.handleError);
  }

  getZernikeDescriptors(
    emdbId: number,
    contourRepresentation: number
  ): Promise<any> {
    return this.httpClient
      .get(
        this.API_URL + "/search/zernike/" + emdbId + "/" + contourRepresentation
      )
      .toPromise()
      .then()
      .catch(this.handleError);
  }

  getZernikeDescriptorsByMapId(
    mapId: number,
    contourRepresentation: number
  ): Promise<any> {
    return this.httpClient
      .get(
        this.API_URL + "/search/zernike/" + mapId + "/" + contourRepresentation
      )
      .toPromise()
      .then()
      .catch(this.handleError);
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
          "/search/results/" +
          emdbId +
          "/" +
          contourRepresentationId +
          "/" +
          isVolumeFilterOn +
          "/" +
          minRes +
          "/" +
          maxRes
      )
      .toPromise()
      .then((data: any) => {
        for (let item of data.results) {
          item.biomolecule.image_url =
            this.API_URL + item.biomolecule.image_url;
        }
        data.path = this.API_URL + data.path;
        return data;
      })
      .catch(this.handleError);
  }

  getSimilarBioMoleculesByMap(
    emdbId: number,
    contourRepresentationId: number,
    isVolumeFilterOn: boolean,
    minRes: number,
    maxRes: number
  ): Promise<any> {
    return this.httpClient
      .get(
        this.API_URL +
          "/search/resultsmap/" +
          emdbId +
          "/" +
          contourRepresentationId +
          "/" +
          isVolumeFilterOn +
          "/" +
          minRes +
          "/" +
          maxRes
      )
      .toPromise()
      .then((data: any) => {
        for (let item of data.results) {
          item.biomolecule.image_url =
            this.API_URL + item.biomolecule.image_url;
        }
        return data;
      })
      .catch(this.handleError);
  }

  getBatchBiomolecules(
    fileId: string,
    contourRepresentationId: number,
    isVolumeFilterOn: boolean,
    topResults: number
  ): Promise<void | BenchmarkResult> {
    return this.httpClient
      .get(
        this.API_URL +
          "/benchmark/query/" +
          fileId +
          "/" +
          contourRepresentationId +
          "/" +
          isVolumeFilterOn +
          "/" +
          topResults
      )
      .toPromise()
      .then((data: BenchmarkResult) => {
        for (let i = 0; i < data.results.length; i++) {
          data.results[i].path = this.API_URL + data.results[i].path;
        }
        data.zipFile = this.API_URL + data.zipFile;
        return data;
      })
      .catch(this.handleError);
  }

  private handleError(error: any) {
    let errMsg = error.message
      ? error.message
      : error.status
      ? `${error.status} - ${error.statusText}`
      : "Server error";
  }
}
