import { TestBed, getTestBed } from "@angular/core/testing";
import {
  HttpClientTestingModule,
  HttpTestingController
} from "@angular/common/http/testing";
import { ErrorHandlerService } from "./error-handler.service";
import { CheckerService } from "./checker.service";

class MockErrorHandlerService {
  handleError(error: any) {}
}

describe("CheckerService", () => {
  let injector: TestBed;
  let service: CheckerService;
  let httpMock: HttpTestingController;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
      providers: [
        CheckerService,
        { provide: ErrorHandlerService, useClass: MockErrorHandlerService }
      ]
    });
    injector = getTestBed();
    service = injector.get(CheckerService);
    httpMock = injector.get(HttpTestingController);
  });

  afterEach(() => {
    httpMock.verify();
  });

  describe("#checkBiomolecule", () => {
    it("should do a get request and return a number", () => {
      const dummyResponse = 1;
      service.checkBiomolecule(5555).then((res: number) => {
        expect(res).toEqual(1);
      });

      const req = httpMock.expectOne(`${service.API_URL}/checker/emdbID/5555`);
      expect(req.request.method).toBe("GET");
      req.flush(dummyResponse);
    });
  });
});
