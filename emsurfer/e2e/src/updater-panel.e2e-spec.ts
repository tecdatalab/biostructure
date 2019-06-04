import { Page } from "./page.po";

describe("Updater and parameters panel tests", () => {
  let page: Page;

  beforeEach(() => {
    page = new Page();
  });

  it("should check if parameters changes are saved", function() {
    page.navigateTo("/admin/parameters");
    let maxVolumeValue = Math.floor(Math.random() * 100 + 1);
    let maxVolume = page.getElementByID("maxVolume");
    maxVolume.click();
    maxVolume.clear();
    maxVolume.sendKeys(maxVolumeValue.toString());
    let hitsInput = page.getElementByID("hitsNumber");
    hitsInput.click();
    hitsInput.clear();
    hitsInput.sendKeys(maxVolumeValue.toString());
    const saveBtn = page.getElementByID("submitButton");
    saveBtn.click();
    page.wait(5000);
    //page.dismissAlert(); //not working in the current version, protractor does not capture alerts
    page.wait(2000);
    maxVolume = page.getElementByID("maxVolume");
    hitsInput = page.getElementByID("hitsNumber");
    expect(maxVolume.getAttribute("value")).toBe(maxVolumeValue.toString());
    expect(hitsInput.getAttribute("value")).toBe(maxVolumeValue.toString());
  });
});
