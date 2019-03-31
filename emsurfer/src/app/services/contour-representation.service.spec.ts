import { TestBed, getTestBed } from "@angular/core/testing";
import { ContourRepresentationService } from "src/app/services/contour-representation.service";
import {
  HttpClientTestingModule,
  HttpTestingController
} from "@angular/common/http/testing";
import { ContourRepresentation } from "../models/contour-representation";
import { ErrorHandlerService } from "./error-handler.service";

class MockErrorHandlerService {
  handleError(error: any) {}
}

describe("ContourRepresentationService", () => {
  let injector: TestBed;
  let service: ContourRepresentationService;
  let httpMock: HttpTestingController;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
      providers: [
        ContourRepresentationService,
        { provide: ErrorHandlerService, useClass: MockErrorHandlerService }
      ]
    });
    injector = getTestBed();
    service = injector.get(ContourRepresentationService);
    httpMock = injector.get(HttpTestingController);
  });

  afterEach(() => {
    httpMock.verify();
  });

  describe("#getContourShapes", () => {
    it("should do a get request and return a ContourRepresentation[]", () => {
      const dummyContourShapes = [
        { id: 0, name: "CR0" },
        { id: 1, name: "CR1" }
      ];
      service.getContourShapes().then((res: ContourRepresentation[]) => {
        expect(res.length).toBe(2);
        expect(res).toEqual(dummyContourShapes);
      });

      const req = httpMock.expectOne(`${service.API_URL}/contour`);
      expect(req.request.method).toBe("GET");
      req.flush(dummyContourShapes);
    });
  });
});
