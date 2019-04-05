import {
  async,
  ComponentFixture,
  TestBed,
  inject
} from "@angular/core/testing";
import { ZernikeDescriptorsListComponent } from "./zernike-descriptors-list.component";

describe("ZernikeDescriptorsListComponent", () => {
  let component: ZernikeDescriptorsListComponent;
  let fixture: ComponentFixture<ZernikeDescriptorsListComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ZernikeDescriptorsListComponent]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(ZernikeDescriptorsListComponent);
    component = fixture.componentInstance;
    component.descriptors = [1, 2, 3];
    component.values = [1.1, 2.2, 3.3];
    fixture.detectChanges();
  });
  it("copy to clipboard button should call cpyToClipboard()", async(() => {
    spyOn(component, "cpyToClipboard");
    const button = fixture.debugElement.nativeElement.querySelector(
      "#cpyToClipboardBtn"
    );
    button.click();
    fixture.whenStable().then(() => {
      expect(component.cpyToClipboard).toHaveBeenCalled();
    });
  }));
});
