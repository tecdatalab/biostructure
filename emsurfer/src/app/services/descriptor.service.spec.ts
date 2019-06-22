import { TestBed, getTestBed } from "@angular/core/testing";
import {
  HttpClientTestingModule,
  HttpTestingController
} from "@angular/common/http/testing";
import { DescriptorService } from "./descriptor.service";
import { ErrorHandlerService } from "./error-handler.service";
import { Descriptor } from "../models/descriptor";
import { DescriptorsList } from "../models/descriptorsList";

class MockErrorHandlerService {
  handleError(error: any) {}
}

describe("DescriptorService", () => {
  let injector: TestBed;
  let service: DescriptorService;
  let httpMock: HttpTestingController;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
      providers: [
        DescriptorService,
        { provide: ErrorHandlerService, useClass: MockErrorHandlerService }
      ]
    });
    injector = getTestBed();
    service = injector.get(DescriptorService);
    httpMock = injector.get(HttpTestingController);
  });

  afterEach(() => {
    httpMock.verify();
  });

  describe("#getDescriptor", () => {
    it("should do a get request and return a Descriptor", () => {
      const dummyDescriptor = new Descriptor();
      dummyDescriptor.emd_entry_id = 1;
      dummyDescriptor.type_descriptor_id = 1;
      dummyDescriptor.numbers = JSON.parse("[1, 2, 3]");
      service.getDescriptor(5555, 1).then((res: Descriptor) => {
        expect(JSON.stringify(res)).toEqual(JSON.stringify(dummyDescriptor));
      });

      const req = httpMock.expectOne(
        `${service.API_URL}/descriptor/zernike/5555/1`
      );
      expect(req.request.method).toBe("GET");
      req.flush(dummyDescriptor);
    });
  });

  describe("#getDescriptorsList", () => {
    it("should do a get request and return a DescriptorList", () => {
      const dummyDescriptorList = new DescriptorsList();
      dummyDescriptorList.descriptors = [];
      dummyDescriptorList.expirationDate = new Date();
      dummyDescriptorList.path = "test";
      service.getDescriptorsList("5555", 1).then((res: DescriptorsList) => {
        expect(JSON.stringify(res)).toEqual(
          JSON.stringify(dummyDescriptorList)
        );
      });

      const req = httpMock.expectOne(
        `${service.API_URL}/descriptor/zernikelist/5555/1`
      );
      expect(req.request.method).toBe("GET");
      req.flush(dummyDescriptorList);
    });
  });
});
