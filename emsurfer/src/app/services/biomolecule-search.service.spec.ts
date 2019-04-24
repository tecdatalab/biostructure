import { TestBed, getTestBed } from "@angular/core/testing";
import {
  HttpClientTestingModule,
  HttpTestingController
} from "@angular/common/http/testing";
import { BiomoleculeSearchService } from "./biomolecule-search.service";
import { Biomolecule } from "../models/biomolecule";
import { SearchResult } from "../models/search-result";
import { BenchmarkResult } from "../models/benchmark-result";
import { ErrorHandlerService } from "./error-handler.service";

class MockErrorHandlerService {
  handleError(error: any) {}
}

describe("BiomoleculeSearchService", () => {
  let injector: TestBed;
  let service: BiomoleculeSearchService;
  let httpMock: HttpTestingController;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
      providers: [
        BiomoleculeSearchService,
        { provide: ErrorHandlerService, useClass: MockErrorHandlerService }
      ]
    });
    injector = getTestBed();
    service = injector.get(BiomoleculeSearchService);
    httpMock = injector.get(HttpTestingController);
  });

  afterEach(() => {
    httpMock.verify();
  });

  describe("#getBiomolecule", () => {
    it("should do a get request and return a Biomolecule", () => {
      const dummyBiomolecule = new Biomolecule();
      dummyBiomolecule.id = 5555;
      service.getBiomolecule(5555).then((res: Biomolecule) => {
        expect(res).toEqual(dummyBiomolecule);
      });
      const req = httpMock.expectOne(`${service.API_URL}/search/5555`);
      expect(req.request.method).toBe("GET");
      req.flush(dummyBiomolecule);
    });
  });

  describe("#getZernikeDescriptors", () => {
    it("should do a get request and return a number[]", () => {
      const dummyDescriptors = [1, 2, 3, 4, 5];
      const emdbId = 5555;
      const contourRepresentation = 0;
      service
        .getZernikeDescriptors(emdbId, contourRepresentation)
        .then((res: number[]) => {
          expect(res).toEqual(dummyDescriptors);
        });
      const req = httpMock.expectOne(
        `${service.API_URL}/search/zernike/5555/0`
      );
      expect(req.request.method).toBe("GET");
      req.flush(dummyDescriptors);
    });
  });

  describe("#getZernikeDescriptorsByMapId", () => {
    it("should do a get request and return a number[]", () => {
      const dummyDescriptors = [1, 2, 3, 4, 5];
      const mapId = 5555;
      const contourRepresentation = 0;
      service
        .getZernikeDescriptorsByMapId(mapId, contourRepresentation)
        .then((res: number[]) => {
          expect(res).toEqual(dummyDescriptors);
        });
      const req = httpMock.expectOne(
        `${service.API_URL}/search/zernike/5555/0`
      );
      expect(req.request.method).toBe("GET");
      req.flush(dummyDescriptors);
    });
  });

  describe("#getSimilarBioMolecules", () => {
    it("should do a get request and return a SearchResult", () => {
      const dummySearchResult = new SearchResult();
      dummySearchResult.results = [];
      dummySearchResult.path = "test";
      const emdbId = 5555;
      const contourRepresentationId = 0;
      const isVolumeFilterOn = true;
      const minRes = "0";
      const maxRes = "0";
      service
        .getSimilarBioMolecules(
          emdbId,
          contourRepresentationId,
          isVolumeFilterOn,
          minRes,
          maxRes
        )
        .then((res: SearchResult) => {
          expect(res).toEqual(dummySearchResult);
        });
      const req = httpMock.expectOne(
        `${
          service.API_URL
        }/search/results/${emdbId}/${contourRepresentationId}/${isVolumeFilterOn}/${minRes}/${maxRes}`
      );
      expect(req.request.method).toBe("GET");
      req.flush(dummySearchResult);
    });
  });

  describe("#getSimilarBioMoleculesByMap", () => {
    it("should do a get request and return a SearchResult", () => {
      const dummySearchResult = new SearchResult();
      dummySearchResult.results = [];
      dummySearchResult.path = "test";
      const contourLevel = 0;
      const contourRepresentationId = 0;
      const filename = "filename";
      const isVolumeFilterOn = true;
      const minRes = "0";
      const maxRes = "0";
      service
        .getSimilarBioMoleculesByMap(
          filename,
          contourRepresentationId,
          contourLevel,
          isVolumeFilterOn,
          minRes,
          maxRes
        )
        .then((res: SearchResult) => {
          expect(res).toEqual(dummySearchResult);
        });
      const req = httpMock.expectOne(
        `${
          service.API_URL
        }/search/resultsmap/${filename}/${contourRepresentationId}/${contourLevel}/${isVolumeFilterOn}/${minRes}/${maxRes}`
      );
      expect(req.request.method).toBe("GET");
      req.flush(dummySearchResult);
    });
  });

  describe("#getBatchBiomolecules", () => {
    it("should do a get request and return a SearchResult", () => {
      const dummyBenchmarkResult = new BenchmarkResult();
      dummyBenchmarkResult.path = "test";
      dummyBenchmarkResult.zipFile = "test";
      dummyBenchmarkResult.results = [];
      const emdbIdList = "5555";
      const contourRepresentationId = 0;
      const isVolumeFilterOn = true;
      const topResults = 20;
      service
        .getBatchBiomolecules(
          emdbIdList,
          contourRepresentationId,
          isVolumeFilterOn,
          topResults
        )
        .then((res: BenchmarkResult) => {
          expect(res).toEqual(dummyBenchmarkResult);
        });
      const req = httpMock.expectOne(
        `${
          service.API_URL
        }/benchmark/query/${emdbIdList}/${contourRepresentationId}/${isVolumeFilterOn}/${topResults}`
      );
      expect(req.request.method).toBe("GET");
      req.flush(dummyBenchmarkResult);
    });
  });
});
