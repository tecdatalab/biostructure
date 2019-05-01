import { AuthServiceConfig } from "angularx-social-login";
import { GoogleLoginProvider } from "angularx-social-login";
export const config = new AuthServiceConfig([
  {
    id: GoogleLoginProvider.PROVIDER_ID,
    provider: new GoogleLoginProvider(
      "930088451397-t93v5v7esll90st0rtc7lvl10uc391cv.apps.googleusercontent.com"
    )
  }
]);
