import { Injectable } from "@angular/core";
import { SearchHistory } from "../models/search-history";
import { HttpClient, HttpHeaders } from "@angular/common/http";
import { UserService } from "../services/user.service";
import { ErrorHandlerService } from "./error-handler.service";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class SearchHistoryService {
  constructor(
    private http: HttpClient,
    private userService: UserService,
    private errorHandlerService: ErrorHandlerService
  ) {}

  readonly API_URL = config.api_url;

  getSearchHistory(): Promise<void | SearchHistory[]> {
    const httpHeaders = new HttpHeaders({
      authorization: this.userService.getStoredAuthToken().token
    });
    return this.http
      .get(this.API_URL + "/history/search", { headers: httpHeaders })
      .toPromise()
      .then((response: SearchHistory[]) => {
        return response;
      })
      .catch(this.errorHandlerService.handleError);
  }
}
