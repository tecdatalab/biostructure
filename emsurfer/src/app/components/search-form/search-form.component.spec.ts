import {
  async,
  ComponentFixture,
  TestBed,
  inject
} from "@angular/core/testing";
import { ReactiveFormsModule } from "@angular/forms";
import { ContourShapeInputComponent } from "../contour-shape-input/contour-shape-input.component";
import { VolumeFilterInputComponent } from "../volume-filter-input/volume-filter-input.component";
import { SearchFormComponent } from "./search-form.component";
import { EmdbIdInputComponent } from "../emdb-id-input/emdb-id-input.component";
import { UploadEmMapComponent } from "../upload-em-map/upload-em-map.component";
import { QueryMethodComponent } from "../query-method/query-method.component";
import { ResolutionFilterComponent } from "../resolution-filter/resolution-filter.component";
import { FileUploadService } from "../../services/file-upload.service";
import { Router } from "@angular/router";
import { HttpClientModule } from "@angular/common/http";
import { CheckerService } from "src/app/services/checker.service";

class MockRouter {
  navigate(urls: string[], extras: string) {
    return true;
  }
}

class MockUploadFilePromise {
  then(func) {
    func(4444);
  }
}

class MockFileUploadService {
  uploadEmMap(urls: string, extras: string) {
    return new MockUploadFilePromise();
  }
}

class MockCheckBiomoleculePromise {
  then(func) {
    return func(5555);
  }
}

class MockCheckerService {
  checkBiomolecule(emdbId) {
    return new MockCheckBiomoleculePromise();
  }
}

describe("SearchFormComponent", () => {
  let component: SearchFormComponent;
  let fixture: ComponentFixture<SearchFormComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      imports: [ReactiveFormsModule, HttpClientModule],
      declarations: [
        SearchFormComponent,
        ContourShapeInputComponent,
        VolumeFilterInputComponent,
        EmdbIdInputComponent,
        UploadEmMapComponent,
        QueryMethodComponent,
        ResolutionFilterComponent
      ],
      providers: [
        { provide: Router, useClass: MockRouter },
        { provide: FileUploadService, useClass: MockFileUploadService },
        { provide: CheckerService, useClass: MockCheckerService }
      ]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(SearchFormComponent);
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

  it("reset function should call this.searchForm.reset(expectedParameter)", async(() => {
    spyOn(component.searchForm, "reset");
    component.reset();
    const expectedParameter = component.defaultFormState;
    const expectedUrl = "result/1884";
    expect(component.searchForm.reset).toHaveBeenCalledWith(expectedParameter);
  }));

  it("submitHandler should call router navigate with expected query params", inject(
    [Router],
    (router: Router) => {
      spyOn(router, "navigate");
      component.submitHandler();
      const expectedQueryParams = {
        contourRepresentation: 1,
        volumeFilter: "On",
        minResolution: null,
        maxResolution: null
      };
      const expectedUrl = "result/1884";
      expect(router.navigate).toHaveBeenCalledWith([expectedUrl], {
        queryParams: expectedQueryParams
      });
    }
  ));

  it("submitHandler should call router navigate with expected query params if the user wants to do a query with a file", inject(
    [Router],
    (router: Router) => {
      spyOn(router, "navigate");
      component.searchForm.patchValue({
        query: {
          search_by_emdb_id: false,
          em_map: {
            filename: "test_file.map"
          }
        }
      });
      component.submitHandler();
      const expectedQueryParams = {
        filename: "test_file.map",
        mapId: 4444,
        contourLevel: 3.14,
        contourRepresentation: 1,
        volumeFilter: "On",
        minResolution: null,
        maxResolution: null
      };
      const expectedUrl = "result/emMap";
      expect(router.navigate).toHaveBeenCalledWith([expectedUrl], {
        queryParams: expectedQueryParams
      });
    }
  ));
});
