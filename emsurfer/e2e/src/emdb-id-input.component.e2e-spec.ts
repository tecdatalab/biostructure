import { EMDBIDPage } from './emdb-id-input.po';

describe("123", () => {
  let page: EMDBIDPage;

  beforeEach(() => {
    page = new EMDBIDPage();
  });

  it("should get the input text", async () => {
    page.navigateTo();

    let inputText = await page.getInput().getAttribute('value');
    expect(inputText).toEqual('1882');
  });
});
