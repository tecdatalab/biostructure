import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import { ErrorHandlerService } from "./error-handler.service";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class CheckerService {
  constructor(
    private httpClient: HttpClient,
    private errorHandlerService: ErrorHandlerService
  ) {}
  readonly API_URL = config.api_url;

  checkBiomolecule(emdbId: number): Promise<void | number> {
    return this.httpClient
      .get(this.API_URL + "/checker/emdbID/" + emdbId)
      .toPromise()
      .then((response: number) => response)
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
