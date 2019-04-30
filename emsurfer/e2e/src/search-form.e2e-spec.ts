import { Page } from "./page.po";

describe("Search form e2e tests", () => {
  let page: Page;

  beforeEach(() => {
    page = new Page();
  });

  it("should prove that EMDB ID input alerts are displayed", function() {
    page.navigateTo("/search");
    const inputText = page.getEmdbIDInput();
    inputText.click();
    inputText.clear();
    let alert = page.getAlert("digits_alert");
    inputText.sendKeys("eeee");
    expect(alert.isDisplayed()).toBe(true);
    inputText.click();
    inputText.clear();
    inputText.sendKeys("222");
    alert = page.getAlert("4_digits_alert");
    expect(alert.isDisplayed()).toBe(true);
  });

  it("should check if submit button is disable", function() {
    page.navigateTo("/search");
    const inputText = page.getEmdbIDInput();
    inputText.click();
    inputText.clear();
    inputText.sendKeys("222");
    expect(page.getButton().isEnabled()).toBe(false);
  });
});
