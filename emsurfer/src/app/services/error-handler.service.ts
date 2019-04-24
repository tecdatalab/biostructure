import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import { Router } from "@angular/router";

@Injectable({
  providedIn: "root"
})
export class ErrorHandlerService {
  constructor(private httpClient: HttpClient, private router: Router) {}

  handleError(error: any) {
    alert("Server response: " + error.error.message + " code: " + error.status);
    if (error.status == 400) {
      this.router.navigate(["/search"]);
    }
  }
}
