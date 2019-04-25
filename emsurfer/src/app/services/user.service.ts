import { Injectable } from "@angular/core";
import config from "../../config.json";
import { HttpClient, HttpHeaders } from "@angular/common/http";
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
    public errorHandlerService: ErrorHandlerService //TODO: service should not be public
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

  changeUserRole(userId: number, role: number) {
    const body = {
      userId: userId,
      role: role
    };
    const httpHeaders = new HttpHeaders({
      authorization: this.getStoredAuthToken().token
    });
    return this.http
      .put(this.API_URL + "/user/admin/changeUserRole", body, {
        headers: httpHeaders
      })
      .toPromise()
      .then((result: any) => {
        return result;
      })
      .catch(this.errorHandlerService.handleError);
  }

  getUserRoles(): Promise<void | UserRole[]> {
    const httpHeaders = new HttpHeaders({
      authorization: this.getStoredAuthToken().token
    });
    return this.http
      .get(this.API_URL + "/user/roles", { headers: httpHeaders })
      .toPromise()
      .then((response: UserRole[]) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  getUsers(): Promise<void | User[]> {
    const httpHeaders = new HttpHeaders({
      authorization: this.getStoredAuthToken().token
    });
    return this.http
      .get(this.API_URL + "/user/users", { headers: httpHeaders })
      .toPromise()
      .then((response: User[]) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }

  checkAdminRole(): Promise<void | boolean> {
    const httpHeaders = new HttpHeaders({
      authorization: this.getStoredAuthToken().token
    });
    return this.http
      .get(this.API_URL + "/user/checkAdminRole", { headers: httpHeaders })
      .toPromise()
      .then((response: boolean) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
