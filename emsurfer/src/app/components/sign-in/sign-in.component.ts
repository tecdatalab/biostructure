import { Component, OnInit } from "@angular/core";
import { AuthService } from "angularx-social-login";
import { GoogleLoginProvider } from "angularx-social-login";
import { SocialUser } from "angularx-social-login";
import { UserService } from "src/app/services/user.service";
import { ErrorHandlerService } from "src/app/services/error-handler.service";
import { Router } from "@angular/router";

@Component({
  selector: "app-sign-in",
  templateUrl: "./sign-in.component.html",
  styleUrls: ["bootstrap-social.css"]
})
export class SignInComponent implements OnInit {
  private loggedIn: boolean;
  private isAdminUser: boolean;
  private errorSubscription;

  constructor(
    private authService: AuthService,
    private userService: UserService,
    private errorHandlerService: ErrorHandlerService,
    private router: Router
  ) {
    this.userService.errorHandlerService = this.errorHandlerService;
    // TO DO: figure out how to share a service between a service and a component in a decent way
  }

  ngOnInit(): void {
    this.errorSubscription = this.errorHandlerService.errorCodeObs.subscribe(
      (errorCode: number) => {
        if (errorCode === 401) {
          this.signOut();
          this.router.navigate(["/home"]);
        }
      }
    );
    this.loggedIn = this.userService.isUserLoggedIn();
    if (this.loggedIn) {
      this.userService.checkAdminRole().then((isAdmin: boolean) => {
        this.isAdminUser = isAdmin;
      });
    }
  }

  ngOnDestroy() {
    this.errorSubscription.unsubscribe();
  }

  signInWithGoogle(): void {
    this.authService
      .signIn(GoogleLoginProvider.PROVIDER_ID)
      .then((user: SocialUser) => {
        this.userService.getAuthToken(user.idToken).then(() => {
          const credential = this.userService.getStoredAuthToken();
          this.loggedIn = credential != null;
          this.userService.checkAdminRole().then((isAdmin: boolean) => {
            this.isAdminUser = isAdmin;
          });
        });
      });
  }

  signOut(): void {
    this.authService.signOut().then(() => {
      this.userService.deleteStoredAuthToken();
      const credential = this.userService.getStoredAuthToken();
      this.loggedIn = credential != null;
    });
  }
}
