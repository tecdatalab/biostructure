import { Injectable } from "@angular/core";
import config from "../../config.json";
import { HttpClient } from "@angular/common/http";
import { ErrorHandlerService } from "./error-handler.service.js";
import { Credential } from "../models/credential.js";
import { UserRole } from "../models/userRole";
import { User } from "../models/user";

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

  isUserLoggedIn() {
    return this.getStoredAuthToken() != null;
  }

  grantAdminRole(userId) {
    const body = {
      token: this.getStoredAuthToken().token,
      userId: userId
    };
    return this.http
      .put(this.API_URL + "/user/admin/grantAdminRole", body)
      .toPromise()
      .then((result: any) => {
        return result;
      })
      .catch(this.errorHandlerService.handleError);
  }

  getUserRoles(): Promise<void | UserRole[]> {
    return this.http
      .get(this.API_URL + "/user/roles")
      .toPromise()
      .then((response: UserRole[]) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  getUsers(): Promise<void | User[]> {
    return this.http
      .get(this.API_URL + "/user/users")
      .toPromise()
      .then((response: User[]) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
