import { Injectable } from "@angular/core";
import { HttpClient } from "@angular/common/http";
import { ErrorHandlerService } from "./error-handler.service";
import config from "../../config.json";
// import * as NGL from "ngl/dist/ngl.esm.js";

@Injectable({
    providedIn: "root"
})

export class BiomoleculeVisualizerService {

    readonly API_URL = config.api_url;

    constructor(
        private httpClient: HttpClient,
        private errorHandlerService: ErrorHandlerService
    ) { }

    //write functions
    getTestStructures(opt: number): Promise<any> {
        return this.httpClient
            .get(this.API_URL + "/visualizer/" + opt)
            .toPromise()
            .then((response: Object) => {
                var structureToShow = []
                Object.keys(response).forEach(item => {
                    structureToShow.push(response[item]);
                });
                return structureToShow;
            })
            .catch(err => {
                this.errorHandlerService.handleError(err);
            });
    }

}
