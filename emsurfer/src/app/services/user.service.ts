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
    console.log(this.getStoredAuthToken());
    let respuesta = this.getStoredAuthToken() != null;
    console.log(respuesta);
    return this.getStoredAuthToken() != null;
  }

  changeUserRole(userId: number, role: number) {
    const body = {
      token: this.getStoredAuthToken().token,
      userId: userId,
      role: role
    };
    return this.http
      .put(this.API_URL + "/user/admin/changeUserRole", body)
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

  checkAdminRole(): Promise<void | boolean> {
    const body = {
      token: this.getStoredAuthToken().token
    };
    return this.http
      .post(this.API_URL + "/user/checkAdminRole", body)
      .toPromise()
      .then((response: boolean) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
