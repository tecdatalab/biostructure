import { browser, by, element } from "protractor";

export class Page {
  navigateTo(path) {
    return browser.get(path);
  }

  navigateToDriver(path) {
    return browser.driver.get(path);
  }

  wait(seconds) {
    browser.driver.sleep(seconds);
    browser.waitForAngular();
  }

  getButton() {
    return element(by.buttonText("Submit"));
  }

  getAlert(idAlert) {
    return element(by.id(idAlert));
  }

  getEmdbIDInput() {
    return element(by.xpath("//input[@id='emdb_id']"));
  }

  getDownloadLink() {
    return element(by.id("download_link"));
  }

  getSignIn() {
    return element(by.id("sign_in"));
  }
}
