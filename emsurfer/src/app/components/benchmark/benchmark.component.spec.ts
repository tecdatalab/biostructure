import {
  async,
  ComponentFixture,
  TestBed,
  inject
} from "@angular/core/testing";
import { ReactiveFormsModule } from "@angular/forms";
import { ContourShapeInputComponent } from "../contour-shape-input/contour-shape-input.component";
import { VolumeFilterInputComponent } from "../volume-filter-input/volume-filter-input.component";
import { BenchmarkComponent } from "./benchmark.component";
import { EmdbIdListComponent } from "../emdb-id-list/emdb-id-list.component";
import { UploadEmdbIdListFileComponent } from "../upload-emdb-id-list-file/upload-emdb-id-list-file.component";
import { Router } from "@angular/router";
import { HttpClientModule } from "@angular/common/http";

class MockRouter {
  navigate(urls: string[], extras: string) {
    return true;
  }
}

describe("BenchmarkComponent", () => {
  let component: BenchmarkComponent;
  let fixture: ComponentFixture<BenchmarkComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      imports: [ReactiveFormsModule, HttpClientModule],
      declarations: [
        ContourShapeInputComponent,
        VolumeFilterInputComponent,
        EmdbIdListComponent,
        UploadEmdbIdListFileComponent,
        BenchmarkComponent
      ],
      providers: [{ provide: Router, useClass: MockRouter }]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(BenchmarkComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it("submit button should call submitHandler()", async(() => {
    spyOn(component, "submitHandler");
    const button = fixture.debugElement.nativeElement.querySelector(
      "#submitButton"
    );
    button.click();
    fixture.whenStable().then(() => {
      expect(component.submitHandler).toHaveBeenCalled();
    });
  }));

  it("reset button should call reset", async(() => {
    spyOn(component, "reset");
    const button = fixture.debugElement.nativeElement.querySelector(
      "#resetButton"
    );
    button.click();
    fixture.whenStable().then(() => {
      expect(component.reset).toHaveBeenCalled();
    });
  }));

  it("submitHandler should call router navigate with expected query params", inject(
    [Router],
    (router: Router) => {
      spyOn(router, "navigate");
      component.submitHandler();
      const expectedQueryParams = {
        emdbIdList: "1010,1884,5502",
        contourRepresentation: 1,
        volumeFilter: "On",
        topResults: 20
      };
      const expectedUrl = "benchmark/results";
      expect(router.navigate).toHaveBeenCalledWith([expectedUrl], {
        queryParams: expectedQueryParams
      });
    }
  ));

  it("EMDB id list Radio buttons should call cbEmdbListChange() method", async(() => {
    spyOn(component, "cbEmdbListChange");
    const rbEmdb = fixture.debugElement.nativeElement.querySelector(
      "#emdb-list"
    );
    rbEmdb.click();
    fixture.whenStable().then(() => {
      expect(component.cbEmdbListChange).toHaveBeenCalled();
    });
  }));

  it("fileupload Radio buttons should call cbEmdbListChange() method", async(() => {
    spyOn(component, "cbEmdbListChange");
    const rbFileupload = fixture.debugElement.nativeElement.querySelector(
      "#fileupload"
    );
    rbFileupload.click();
    fixture.whenStable().then(() => {
      expect(component.cbEmdbListChange).toHaveBeenCalled();
    });
  }));
});
