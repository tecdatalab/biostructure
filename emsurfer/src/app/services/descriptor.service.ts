import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import { Descriptor } from "../models/descriptor";
import { DescriptorsList } from "../models/descriptorsList";
import { ErrorHandlerService } from "./error-handler.service";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class DescriptorService {
  constructor(
    private httpClient: HttpClient,
    private errorHandlerService: ErrorHandlerService
  ) {}

  readonly API_URL = config.api_url;

  getDescriptor(
    emdbId: number,
    contourRepresentation: number
  ): Promise<void | Descriptor> {
    return this.httpClient
      .get(
        this.API_URL +
          "/descriptor/zernike/" +
          emdbId +
          "/" +
          contourRepresentation
      )
      .toPromise()
      .then((response: Descriptor) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  getDescriptorsList(
    emdbList: string,
    contourRepresentation: number
  ): Promise<void | DescriptorsList> {
    return this.httpClient
      .get(
        this.API_URL +
          "/descriptor/zernikelist/" +
          emdbList +
          "/" +
          contourRepresentation
      )
      .toPromise()
      .then((response: DescriptorsList) => {
        response.path = this.API_URL + response.path;
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
