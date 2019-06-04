import { Injectable } from "@angular/core";
import { HttpClient, HttpHeaders } from "@angular/common/http";
import { Update } from "../models/update";
import { UserService } from "../services/user.service";
import { ErrorHandlerService } from "./error-handler.service";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class UpdateService {
  constructor(
    private http: HttpClient,
    private userService: UserService,
    private errorHandlerService: ErrorHandlerService
  ) {}

  readonly API_URL = config.api_url;

  getLastUpdate(): Promise<void | Update> {
    const httpHeaders = new HttpHeaders({
      authorization: this.userService.getStoredAuthToken().token
    });
    return this.http
      .get(this.API_URL + "/update", { headers: httpHeaders })
      .toPromise()
      .then((response: Update) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  forceUpdate() {
    const httpHeaders = new HttpHeaders({
      authorization: this.userService.getStoredAuthToken().token
    });
    return this.http
      .get(this.API_URL + "/update/forceUpdater", { headers: httpHeaders })
      .toPromise()
      .then((response: Update) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
