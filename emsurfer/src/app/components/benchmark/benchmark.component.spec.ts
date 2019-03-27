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

describe("SearchFormComponent", () => {
  let component: SearchFormComponent;
  let fixture: ComponentFixture<SearchFormComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      imports: [ReactiveFormsModule, RouterTestingModule],
      declarations: [
        SearchFormComponent,
        HeaderComponent,
        ContourShapeInputComponent,
        VolumeFilterInputComponent,
        EmdbIdInputComponent,
        UploadEmMapComponent,
        QueryMethodComponent,
        ResolutionFilterComponent
      ],
      providers: [
        { provide: Router, useClass: MockRouter },
        { provide: FileUploadService, useClass: MockFileUploadService }
      ]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(SearchFormComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });
});
