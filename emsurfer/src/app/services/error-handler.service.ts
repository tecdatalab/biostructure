import { Injectable } from "@angular/core";
import { Router } from "@angular/router";
import { BehaviorSubject } from "rxjs";

@Injectable({
  providedIn: "root"
})
export class ErrorHandlerService {
  private errorCode = new BehaviorSubject<number>(null);
  errorCodeObs = this.errorCode.asObservable();
  constructor(private router: Router) {}

  handleError(error: any) {
    if (error.status === 401) {
      alert("Your session has expired, please log in and try again");
    } else {
      alert(
        "Server response: " + error.error.message + " code: " + error.status
      );
    }
    if (error.status === 400) {
      this.router.navigate(["/search"]);
    }
    this.errorCode.next(error.status);
  }
}
