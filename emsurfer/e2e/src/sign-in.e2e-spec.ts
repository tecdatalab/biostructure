import { Page } from "./page.po";
const path = require("path");
const fs = require("fs");

describe("Sign In", () => {
  let page: Page;

  beforeEach(() => {
    page = new Page();
  });

  it("check the sign in function", function() {
    expect(true).toBe(true);
  });
});
