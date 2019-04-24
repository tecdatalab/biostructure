import { TestBed, getTestBed } from "@angular/core/testing";
import {
  HttpClientTestingModule,
  HttpTestingController
} from "@angular/common/http/testing";
import { ErrorHandlerService } from "./error-handler.service";
import { Credential } from "../models/credential";
import { User } from "../models/user";
import { UserService } from "./user.service";
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

fdescribe("UserService", () => {
  let injector: TestBed;
  let service: UserService;
  let httpMock: HttpTestingController;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
      providers: [
        UserService,
        { provide: Router, useClass: MockRouter },
        { provide: ErrorHandlerService, useClass: MockErrorHandlerService }
      ]
    });
    injector = getTestBed();
    service = injector.get(UserService);
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

  describe("#getAuthToken", () => {
    it("should do a post request and should store credential in local storage", () => {
      const googleTokenId = "googleToken";
      const expectedBody = {
        tokenId: googleTokenId
      };
      service.getAuthToken(googleTokenId).then(() => {
        const storedCredential = window.localStorage.getItem("credential");
        expect(storedCredential).toEqual(JSON.stringify(dummyCredential));
      });
      const req = httpMock.expectOne(`${service.API_URL}/user/auth/token`);
      expect(req.request.method).toBe("POST");
      expect(req.request.body).toEqual(expectedBody);
      req.flush(dummyCredential);
    });
  });

  describe("#getStoredAuthToken", () => {
    it("should return expected token stored in localstorage", () => {
      window.localStorage.setItem(
        "credential",
        JSON.stringify(dummyCredential)
      );
      const actualCredential = service.getStoredAuthToken();
      expect(JSON.stringify(actualCredential)).toEqual(
        JSON.stringify(dummyCredential)
      );
      expect(window.localStorage.getItem).toHaveBeenCalledWith("credential");
    });
  });

  describe("#deleteStoredAuthToken", () => {
    it("should credential stored in localstorage", () => {
      service.deleteStoredAuthToken();
      expect(window.localStorage.getItem("credential")).toEqual(null);
      expect(window.localStorage.removeItem).toHaveBeenCalledWith("credential");
    });
  });

  describe("#isUserLoggedIn", () => {
    it("should return true if the user is logged in", () => {
      window.localStorage.setItem(
        "credential",
        JSON.stringify(dummyCredential)
      );
      expect(service.isUserLoggedIn()).toEqual(true);
    });
    it("should return false if the user is not logged in", () => {
      expect(service.isUserLoggedIn()).toEqual(false);
    });
  });

  describe("#changeUserRole", () => {
    it("should do a put request and should return the updated user", () => {
      window.localStorage.setItem(
        "credential",
        JSON.stringify(dummyCredential)
      );
      const expectedHeader = dummyCredential.token;
      const expectedUser = {
        id: dummyCredential.user.id,
        role: 2
      };

      service.changeUserRole(+dummyCredential.user.id, 2).then(result => {
        expect(JSON.stringify(result)).toEqual(JSON.stringify(expectedUser));
      });
      const req = httpMock.expectOne(
        `${service.API_URL}/user/admin/changeUserRole`
      );
      expect(req.request.method).toBe("PUT");
      expect(req.request.headers.get("authorization")).toEqual(expectedHeader);
      req.flush(expectedUser);
    });
  });

  describe("#getUserRoles", () => {
    it("should do a get request and should return a list of user roles", () => {
      window.localStorage.setItem(
        "credential",
        JSON.stringify(dummyCredential)
      );
      const expectedHeader = dummyCredential.token;
      const expectedRoles = [
        { id: 1, role: "role1" },
        {
          id: 2,
          role: "role2"
        }
      ];
      service.getUserRoles().then(result => {
        expect(JSON.stringify(result)).toEqual(JSON.stringify(expectedRoles));
      });
      const req = httpMock.expectOne(`${service.API_URL}/user/roles`);
      expect(req.request.method).toBe("GET");
      expect(req.request.headers.get("authorization")).toEqual(expectedHeader);
      req.flush(expectedRoles);
    });
  });

  describe("#getUsers", () => {
    it("should do a get request and should return a list of users", () => {
      window.localStorage.setItem(
        "credential",
        JSON.stringify(dummyCredential)
      );
      const expectedHeader = dummyCredential.token;
      const expectedUsers = [
        {
          id: "1",
          name: "name1",
          email: "email1",
          role: 1
        },
        {
          id: "2",
          name: "name2",
          email: "email2",
          role: 2
        }
      ];
      service.getUsers().then(results => {
        expect(JSON.stringify(results)).toEqual(JSON.stringify(expectedUsers));
      });
      const req = httpMock.expectOne(`${service.API_URL}/user/users`);
      expect(req.request.method).toBe("GET");
      expect(req.request.headers.get("authorization")).toEqual(expectedHeader);
      req.flush(expectedUsers);
    });
  });

  describe("#checkAdminRole", () => {
    it("should do a get request and should return a boolean", () => {
      window.localStorage.setItem(
        "credential",
        JSON.stringify(dummyCredential)
      );
      const expectedHeader = dummyCredential.token;
      const expectedResult = true;
      service.checkAdminRole().then(result => {
        expect(!!result).toEqual(expectedResult);
      });
      const req = httpMock.expectOne(`${service.API_URL}/user/checkAdminRole`);
      expect(req.request.method).toBe("GET");
      expect(req.request.headers.get("authorization")).toEqual(expectedHeader);
      req.flush(+expectedResult);
    });
  });
});
