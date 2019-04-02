import { browser, by, element } from "protractor";

export class EMDBIDPage {
  navigateTo() {
    return browser.get("/");
  }

  getInput() {
    return element(
      by.css(
        '[placeholder = "EMDB entry ID"]'
      )
    );
  }
}
