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

  getDescriptor(emdbId: number): Promise<void | Descriptor> {
    return this.httpClient
      .get(this.API_URL + "/descriptor/zernike/" + emdbId)
      .toPromise()
      .then((response: Descriptor) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  getDescriptorsList(emdbList: string): Promise<void | DescriptorsList> {
    return this.httpClient
      .get(this.API_URL + "/descriptor/zernikelist/" + emdbList)
      .toPromise()
      .then((response: DescriptorsList) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
