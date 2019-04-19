import { Component } from "@angular/core";
import { AuthService } from "angularx-social-login";
import { GoogleLoginProvider } from "angularx-social-login";
import { SocialUser } from "angularx-social-login";
import { UserService } from "src/app/services/user.service";

@Component({
  selector: "app-sign-in",
  templateUrl: "./sign-in.component.html",
  styleUrls: ["bootstrap-social.css"]
})
export class SignInComponent {
  private loggedIn: boolean;

  constructor(
    private authService: AuthService,
    private userService: UserService
  ) {}

  signInWithGoogle(): void {
    this.authService
      .signIn(GoogleLoginProvider.PROVIDER_ID)
      .then((user: SocialUser) => {
        this.userService.getAuthToken(user.idToken).then(() => {
          const credential = this.userService.getStoredAuthToken();
          this.loggedIn = credential != null;
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
