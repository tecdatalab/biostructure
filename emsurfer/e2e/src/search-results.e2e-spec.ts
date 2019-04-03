import { Page } from './page.po';
const path = require('path');
const downloadsPath = path.resolve(__dirname, '../downloads');
const fs = require('fs');

describe('Download test', () => {
  let page: Page;

  beforeEach(() => {
    page = new Page();
  });

  it('check if the result file was downloaded', function() {
    page.navigateTo('/result/1884?contourRepresentation=0&volumeFilter=On&minResolution=0&maxResolution=12');
    page.wait(4000);
    const downloadLink =  page.getDownloadLink();
    downloadLink.click();
    console.log(downloadsPath);
    page.wait(4000);
    const fileExists = fs.existsSync(downloadsPath);
    expect(fileExists).toBe(true);
  });

});
