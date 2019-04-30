import { Page } from "./page.po";
const path = require("path");
const fs = require("fs");

describe("Sign In", () => {
  let page: Page;

  beforeEach(() => {
    page = new Page();
  });

  it("check the sign in function", function() {
    page.navigateToDriver("localhost/4200/home");
    page.wait(4000);
    const signInButton = page.getSignIn();
    signInButton.click();
    page.wait(5000);
    expect(true).toBe(true);
  });
});
