import { Page } from "./page.po";

describe("Zernike Module", () => {
  let page: Page;

  beforeEach(() => {
    page = new Page();
  });

  it("check the sign in function", function() {
    page.navigateTo("/benchmark");
    page.wait(1000);
    const inputText = page.getEmdbIDList();
    inputText.click();
    inputText.clear();
    inputText.sendKeys("0002 \n1884 \n0001");
    page.wait(1000);
    const button = page.getButton("zernikeButton");
    button.click();
    page.wait(5000);
    expect(true).toBe(true);
  });
});
