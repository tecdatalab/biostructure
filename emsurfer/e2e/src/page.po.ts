import { browser, by, element, protractor } from "protractor";

export class Page {
  dismissAlert() {
    browser
      .switchTo()
      .alert()
      .accept();
  }

  navigateTo(path) {
    return browser.get(path);
  }

  getCurrentUrl() {
    return browser.getCurrentUrl();
  }

  navigateToDriver(path) {
    return browser.driver.get(path);
  }

  wait(seconds) {
    browser.driver.sleep(seconds);
    browser.waitForAngular();
  }

  getSubmitButton() {
    return element(by.buttonText("Submit"));
  }

  getButton(buttonId) {
    return element(by.id(buttonId));
  }

  getAlert(idAlert) {
    return element(by.id(idAlert));
  }

  getDownloadZernike() {
    return element(by.linkText("Download text results here"));
  }

  getEmdbIDInput() {
    return element(by.xpath("//input[@id='emdb_id']"));
  }

  getEmdbIDList() {
    return element(by.id("emdbIdList"));
  }

  getDownloadLink() {
    return element(by.id("download_link"));
  }

  getSignIn() {
    return element(by.id("sign_in"));
  }

  getElementByID(elementID) {
    return element(by.id(elementID));
  }
}
