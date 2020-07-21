import { TestBed, getTestBed } from "@angular/core/testing";
import {
    HttpClientTestingModule,
    HttpTestingController
} from "@angular/common/http/testing";
import { BiomoleculeVisualizerService } from "./biomolecule-visualizer.service";
import { ErrorHandlerService } from "./error-handler.service";

class MockErrorHandlerService {
    handleError(error: any) { }
}

describe("BiomoleculeVisualizerService", () => {

    let injector: TestBed;
    let service: BiomoleculeVisualizerService;
    let httpMock: HttpTestingController;

    
})

