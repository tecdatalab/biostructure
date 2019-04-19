import { Injectable } from "@angular/core";
import config from "../../config.json";
import { HttpClient } from "@angular/common/http";
import { ErrorHandlerService } from "./error-handler.service.js";
import { Credential } from "../models/credential.js";

@Injectable({
  providedIn: "root"
})
export class UserService {
  API_URL = config.api_url;
  constructor(
    private http: HttpClient,
    private errorHandlerService: ErrorHandlerService
  ) {}

  getAuthToken(googleIdToken: string): Promise<void | void> {
    const body = {
      tokenId: googleIdToken
    };
    return this.http
      .post(this.API_URL + "/user/auth/token", body)
      .toPromise()
      .then((credential: Credential) => {
        localStorage.setItem("credential", JSON.stringify(credential));
      })
      .catch(this.errorHandlerService.handleError);
  }

  getStoredAuthToken(): Credential {
    const storedCredential = JSON.parse(localStorage.getItem("credential"));
    return storedCredential;
  }

  deleteStoredAuthToken(): void {
    localStorage.removeItem("credential");
  }
}
