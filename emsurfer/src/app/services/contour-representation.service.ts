import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import { ContourRepresentation } from "../models/contour-representation";
import { ErrorHandlerService } from "./error-handler.service";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class ContourRepresentationService {
  constructor(
    private httpClient: HttpClient,
    private errorHandlerService: ErrorHandlerService
  ) {}

  readonly API_URL = config.api_url;

  getContourShapes(): Promise<void | ContourRepresentation[]> {
    return this.httpClient
      .get(this.API_URL + "/contour")
      .toPromise()
      .then((response: ContourRepresentation[]) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
