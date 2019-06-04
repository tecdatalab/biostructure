import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import { ErrorHandlerService } from "./error-handler.service";
import config from "../../config.json";
import { Statistic } from "../models/statistic";

@Injectable({
  providedIn: "root"
})
export class StatisticsService {
  API_URL = config.api_url;

  constructor(
    private httpClient: HttpClient,
    private errorHandlerService: ErrorHandlerService
  ) {}

  getStatistics(): Promise<void | Statistic[]> {
    return this.httpClient
      .get(this.API_URL + "/statistics/")
      .toPromise()
      .then((response: Statistic[]) => {
        return response;
      })
      .catch(err => {
        this.errorHandlerService.handleError(err);
      });
  }
}
