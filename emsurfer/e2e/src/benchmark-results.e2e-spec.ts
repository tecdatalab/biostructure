import { Page } from './page.po';
const path = require('path');
const downloadsPath = path.resolve(__dirname, '../downloads/results.zip');
const fs = require('fs');

describe('Download benchmark test', () => {
  let page: Page;

  beforeEach(() => {
    page = new Page();
  });

  it('check if the benchmark result file was downloaded and the file format (zip)', function() {
    page.navigateTo('/benchmark/results?emdbIdList=1010,1884,5502&contourRepresentation=0&volumeFilter=On&topResults=20');
    page.wait(4000);
    const downloadLink =  page.getDownloadLink();
    downloadLink.click();
    page.wait(3000);
    const fileExists = fs.existsSync(downloadsPath);
    page.wait(2000);
    expect(fileExists).toBe(true);
  });

});
