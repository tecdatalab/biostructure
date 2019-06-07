import { Page } from "./page.po";

describe("Search history test", () => {
  let page: Page;

  beforeEach(() => {
    page = new Page();
  });

  it("check if new search is registered", function() {
    page.navigateTo("/reference");
    page.wait(15000);
    page.navigateTo("/search");
    const inputText = page.getEmdbIDInput();
    inputText.click();
    inputText.clear();
    inputText.sendKeys("0008");
    page.wait(500);
    page.getButton("submitButton").click();
    page.wait(2000);
    page.navigateTo("/history");
    const inputDate = page.getElementByID("minDate");
    let newDate = new Date();
    inputDate.sendKeys(
      ("0" + newDate.getDate()).slice(-2) +
        ("0" + newDate.getMonth()).slice(-2) +
        newDate.getFullYear()
    );
    page.wait(1000);
    const searchBtn = page.getButton("searchBtn");
    page.wait(5000);
    expect(searchBtn.isEnabled()).toBe(true);
  });

  it("check if search button works", function() {
    page.navigateTo("/history");
    page.wait(5000);
    const searchBtn = page.getButton("searchBtn");
    searchBtn.click();
    page.wait(5000);
    expect(true).toBe(true);
    expect(page.getCurrentUrl()).toContain("http://localhost:4200/result");
  });
});
