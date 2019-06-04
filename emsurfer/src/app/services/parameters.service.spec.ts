import { TestBed, getTestBed } from "@angular/core/testing";
import {
  HttpClientTestingModule,
  HttpTestingController
} from "@angular/common/http/testing";
import { ParametersService } from "./parameters.service";
import { ErrorHandlerService } from "./error-handler.service";
import { Parameters } from "../models/parameters";
import { Credential } from "../models/credential";
import { User } from "../models/user";
import { Router } from "@angular/router";

class MockErrorHandlerService {
  handleError(error: any) {}
}
class MockRouter {
  navigate(urls: string[], extras: string) {
    return true;
  }
}

const dummyCredential = new Credential();
dummyCredential.token = "token";
dummyCredential.user = new User();
dummyCredential.user.id = "id";

describe("ParametersService", () => {
  let injector: TestBed;
  let service: ParametersService;
  let httpMock: HttpTestingController;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
      providers: [
        ParametersService,
        { provide: Router, useClass: MockRouter },
        { provide: ErrorHandlerService, useClass: MockErrorHandlerService }
      ]
    });
    injector = getTestBed();
    service = injector.get(ParametersService);
    httpMock = injector.get(HttpTestingController);

    let store = {};

    spyOn(window.localStorage, "getItem").and.callFake(
      (key: string): String => {
        return store[key] || null;
      }
    );
    spyOn(window.localStorage, "removeItem").and.callFake(
      (key: string): void => {
        delete store[key];
      }
    );
    spyOn(window.localStorage, "setItem").and.callFake(
      (key: string, value: string): string => {
        return (store[key] = <string>value);
      }
    );
  });

  afterEach(() => {
    httpMock.verify();
  });

  describe("#getParameters", () => {
    it("should do a get request and return a parameters object", () => {
      const dummyParameters = new Parameters();
      dummyParameters.hits_number = 20;
      dummyParameters.update_rate = 20;
      dummyParameters.volume_filter_max = 0.5;
      dummyParameters.volume_filter_min = 0.5;
      window.localStorage.setItem(
        "credential",
        JSON.stringify(dummyCredential)
      );
      const expectedHeader = dummyCredential.token;
      service.getParameters().then((res: Parameters) => {
        expect(JSON.stringify(res)).toEqual(JSON.stringify(dummyParameters));
      });

      const req = httpMock.expectOne(`${service.API_URL}/parameters`);
      expect(req.request.method).toBe("GET");
      expect(req.request.headers.get("authorization")).toEqual(expectedHeader);
      req.flush(dummyParameters);
    });
  });

  describe("#setParameters", () => {
    it("should do a get request", () => {
      window.localStorage.setItem(
        "credential",
        JSON.stringify(dummyCredential)
      );
      const expectedHeader = dummyCredential.token;
      service.setParameters(0.5, 0.6, 5, 6).then((res: any) => {});

      const req = httpMock.expectOne(
        `${service.API_URL}/parameters/set/0.5/0.6/5/6`
      );
      expect(req.request.method).toBe("GET");
      expect(req.request.headers.get("authorization")).toEqual(expectedHeader);
      req.flush(null);
    });
  });
});
