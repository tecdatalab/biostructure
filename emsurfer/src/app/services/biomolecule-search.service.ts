import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import { Biomolecule } from "../models/biomolecule";
import { BenchmarkResult } from "../models/benchmark-result";
import { SearchResult } from "../models/search-result";
import { ErrorHandlerService } from "./error-handler.service";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class BiomoleculeSearchService {
  constructor(
    private httpClient: HttpClient,
    private errorHandlerService: ErrorHandlerService
  ) {}

  readonly API_URL = config.api_url;

  getBiomolecule(emdbId: number): Promise<void | Biomolecule> {
    return this.httpClient
      .get(this.API_URL + "/search/" + emdbId)
      .toPromise()
      .then((response: Biomolecule) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  getZernikeDescriptors(
    emdbId: number,
    contourRepresentation: number
  ): Promise<void | number[]> {
    return this.httpClient
      .get(
        this.API_URL + "/search/zernike/" + emdbId + "/" + contourRepresentation
      )
      .toPromise()
      .then((response: number[]) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  getZernikeDescriptorsByMapId(
    mapId: number,
    contourRepresentation: number
  ): Promise<void | number[]> {
    return this.httpClient
      .get(
        this.API_URL + "/search/zernike/" + mapId + "/" + contourRepresentation
      )
      .toPromise()
      .then((response: number[]) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  getSimilarBioMolecules(
    emdbId: number,
    contourRepresentationId: number,
    isVolumeFilterOn: boolean,
    minRes: string,
    maxRes: string
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
      .then((data: SearchResult) => {
        data.path = this.API_URL + data.path;
        return data;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  getSimilarBioMoleculesByMap(
    filename: string,
    contourRepresentationId: number,
    contourLevel: number,
    isVolumeFilterOn: boolean,
    minRes: string,
    maxRes: string
  ): Promise<any> {
    return this.httpClient
      .get(
        this.API_URL +
          "/search/resultsmap/" +
          filename +
          "/" +
          contourRepresentationId +
          "/" +
          contourLevel +
          "/" +
          isVolumeFilterOn +
          "/" +
          minRes +
          "/" +
          maxRes
      )
      .toPromise()
      .then((data: SearchResult) => {
        data.path = this.API_URL + data.path;
        return data;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  getBatchBiomolecules(
    emdbIdList: string,
    contourRepresentationId: number,
    isVolumeFilterOn: boolean,
    topResults: number
  ): Promise<void | BenchmarkResult> {
    return this.httpClient
      .get(
        this.API_URL +
          "/benchmark/query/" +
          emdbIdList +
          "/" +
          contourRepresentationId +
          "/" +
          isVolumeFilterOn +
          "/" +
          topResults
      )
      .toPromise()
      .then((data: BenchmarkResult) => {
        data.zipFile = this.API_URL + data.zipFile;
        return data;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
