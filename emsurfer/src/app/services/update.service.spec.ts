import { TestBed, getTestBed } from "@angular/core/testing";
import {
  HttpClientTestingModule,
  HttpTestingController
} from "@angular/common/http/testing";
import { UpdateService } from "./update.service";
import { ErrorHandlerService } from "./error-handler.service";
import { Credential } from "../models/credential";
import { User } from "../models/user";
import { Router } from "@angular/router";
import { Update } from "../models/update";

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

describe("UpdateService", () => {
  let injector: TestBed;
  let service: UpdateService;
  let httpMock: HttpTestingController;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
      providers: [
        UpdateService,
        { provide: Router, useClass: MockRouter },
        { provide: ErrorHandlerService, useClass: MockErrorHandlerService }
      ]
    });
    injector = getTestBed();
    service = injector.get(UpdateService);
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

  describe("#getLastUpdate", () => {
    it("should do a get request and return an update object", () => {
      const dummyUpdate = new Update();
      dummyUpdate.last_update = new Date("00-00-00");
      window.localStorage.setItem(
        "credential",
        JSON.stringify(dummyCredential)
      );
      const expectedHeader = dummyCredential.token;
      service.getLastUpdate().then((res: Update) => {
        expect(JSON.stringify(res)).toEqual(JSON.stringify(dummyUpdate));
      });

      const req = httpMock.expectOne(`${service.API_URL}/update`);
      expect(req.request.method).toBe("GET");
      expect(req.request.headers.get("authorization")).toEqual(expectedHeader);
      req.flush(dummyUpdate);
    });
  });

  describe("#forceUpdate", () => {
    it("should do a get request and return an update object", () => {
      const dummyUpdate = new Update();
      dummyUpdate.last_update = new Date("00-00-00");
      window.localStorage.setItem(
        "credential",
        JSON.stringify(dummyCredential)
      );
      const expectedHeader = dummyCredential.token;
      service.forceUpdate().then((res: Update) => {
        expect(JSON.stringify(res)).toEqual(JSON.stringify(dummyUpdate));
      });

      const req = httpMock.expectOne(`${service.API_URL}/update/forceUpdater`);
      expect(req.request.method).toBe("GET");
      expect(req.request.headers.get("authorization")).toEqual(expectedHeader);
      req.flush(dummyUpdate);
    });
  });
});
