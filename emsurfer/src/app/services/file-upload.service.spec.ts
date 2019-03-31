import { TestBed, getTestBed } from "@angular/core/testing";
import { ContourRepresentationService } from "src/app/services/contour-representation.service";
import {
  HttpClientTestingModule,
  HttpTestingController
} from "@angular/common/http/testing";
import { FileUploadService } from "./file-upload.service";
import { ErrorHandlerService } from "./error-handler.service";

class MockErrorHandlerService {
  handleError(error: any) {}
}

describe("FileUploadService", () => {
  let injector: TestBed;
  let service: FileUploadService;
  let httpMock: HttpTestingController;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
      providers: [
        FileUploadService,
        { provide: ErrorHandlerService, useClass: MockErrorHandlerService }
      ]
    });
    injector = getTestBed();
    service = injector.get(FileUploadService);
    httpMock = injector.get(HttpTestingController);
  });

  afterEach(() => {
    httpMock.verify();
  });

  describe("#uploadEmMap", () => {
    it("should do a post request and return a number", () => {
      const dummyId = 5;
      const fileBase64 = "base64";
      const filename = "filename";
      const expectedBody = {
        file: fileBase64,
        filename: filename
      };
      service.uploadEmMap("base64", "filename").then((res: number) => {
        expect(res).toEqual(dummyId);
      });

      const req = httpMock.expectOne(`${service.API_URL}/upload/EmMap`);
      expect(req.request.method).toBe("POST");
      expect(req.request.body).toEqual(expectedBody);
      req.flush(dummyId);
    });
  });
});
