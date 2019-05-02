import { Page } from "./page.po";
const path = require("path");
const downloadsPath = path.resolve(__dirname, "../downloads/result.hit");
const fs = require("fs");

describe("Download test", () => {
  let page: Page;

  beforeEach(() => {
    page = new Page();
  });

  it("check if the search result file was downloaded and the file format (hit)", function() {
    page.navigateTo(
      "/result/0002?contourRepresentation=0&volumeFilter=On&minResolution=0&maxResolution=12"
    );
    page.wait(4000);
    const downloadLink = page.getDownloadLink();
    downloadLink.click();
    page.wait(5000);
    const fileExists = fs.existsSync(downloadsPath);
    page.wait(5000);
    expect(fileExists).toBe(true);
  });
});
