import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import { ErrorHandlerService } from "./error-handler.service";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class FileUploadService {
  constructor(
    private http: HttpClient,
    private errorHandlerService: ErrorHandlerService
  ) {}

  API_URL = config.api_url;

  uploadEmMap(base64File: string, file_name: string): Promise<void | number> {
    const body = {
      file: base64File,
      filename: file_name
    };
    return this.http
      .post(this.API_URL + "/upload/EmMap", body)
      .toPromise()
      .then((data: number) => data)
      .catch(this.errorHandlerService.handleError);
  }
}
