import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import { ContourRepresentation } from "../models/contour-representation";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class ContourRepresentationService {
  constructor(private httpClient: HttpClient) {}

  readonly API_URL = config.api_url;

  getContourShapes(): Promise<void | ContourRepresentation[]> {
    return this.httpClient
      .get(this.API_URL + "/contour")
      .toPromise()
      .then((response: ContourRepresentation[]) => {
        console.log(response);
        return response;
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
