import { Injectable } from "@angular/core";
import { HttpClient, HttpHeaders } from "@angular/common/http";
import { Parameters } from "../models/parameters";
import { UserService } from "../services/user.service";
import { ErrorHandlerService } from "./error-handler.service";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class ParametersService {
  constructor(
    private http: HttpClient,
    private userService: UserService,
    private errorHandlerService: ErrorHandlerService
  ) {}

  readonly API_URL = config.api_url;

  getParameters(): Promise<void | Parameters> {
    const httpHeaders = new HttpHeaders({
      authorization: this.userService.getStoredAuthToken().token
    });
    return this.http
      .get(this.API_URL + "/parameters", { headers: httpHeaders })
      .toPromise()
      .then((response: Parameters) => {
        return response;
      });
  }

  setParameters(
    volumeMin: number,
    volumeMax: number,
    hits: number,
    update: number
  ) {
    const httpHeaders = new HttpHeaders({
      authorization: this.userService.getStoredAuthToken().token
    });
    this.http
      .get(
        this.API_URL +
          "/parameters/set/" +
          volumeMin +
          "/" +
          volumeMax +
          "/" +
          hits +
          "/" +
          update,
        { headers: httpHeaders }
      )
      .toPromise()
      .then()
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
