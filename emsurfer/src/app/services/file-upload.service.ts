import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import config from "../../config.json";

@Injectable({
  providedIn: "root"
})
export class FileUploadService {
  constructor(private http: HttpClient) {}

  API_URL = config.api_url;

  uploadIdListFile(base64File: string) {
    return 5555;
  }

  uploadIdListString(idList: string) {
    return 1234;
  }

  uploadEmMap(base64File: string): Promise<void | number> {
    console.log("amigos")
    return this.http
      .post(this.API_URL + "/upload/EmMap", base64File)
      .toPromise()
      .then((data: number) => data)
      .catch();
  }
}
