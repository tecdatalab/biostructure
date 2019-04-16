import { BrowserModule } from "@angular/platform-browser";
import { NgModule, ErrorHandler } from "@angular/core";
import { ReactiveFormsModule } from "@angular/forms";
import { HttpClientModule } from "@angular/common/http";

import { routes } from "./app.router";

import { AppComponent } from "./app.component";
import { HeaderComponent } from "./components/header/header.component";
import { ContourShapeInputComponent } from "./components/contour-shape-input/contour-shape-input.component";
import { VolumeFilterInputComponent } from "./components/volume-filter-input/volume-filter-input.component";
import { SearchFormComponent } from "./components/search-form/search-form.component";
import { EmdbIdInputComponent } from "./components/emdb-id-input/emdb-id-input.component";
import { UploadEmMapComponent } from "./components/upload-em-map/upload-em-map.component";
import { QueryMethodComponent } from "./components/query-method/query-method.component";
import { ResolutionFilterComponent } from "./components/resolution-filter/resolution-filter.component";
import { BiomoleculeComponent } from "./components/biomolecule/biomolecule.component";
import { BiomoleculesTableComponent } from "./components/biomolecules-table/biomolecules-table.component";
import { SearchResultComponent } from "./components/search-result/search-result.component";
import { ZernikeDescriptorsListComponent } from "./components/zernike-descriptors-list/zernike-descriptors-list.component";
import { BenchmarkComponent } from "./components/benchmark/benchmark.component";
import { EmdbIdListComponent } from "./components/emdb-id-list/emdb-id-list.component";
import { UploadEmdbIdListFileComponent } from "./components/upload-emdb-id-list-file/upload-emdb-id-list-file.component";
import { BenchmarkResultsComponent } from "./components/benchmark-results/benchmark-results.component";
import { HomeComponent } from './components/home/home.component';
import { ReferenceComponent } from './components/reference/reference.component';

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    ContourShapeInputComponent,
    VolumeFilterInputComponent,
    SearchFormComponent,
    EmdbIdInputComponent,
    UploadEmMapComponent,
    QueryMethodComponent,
    ResolutionFilterComponent,
    BiomoleculeComponent,
    BiomoleculesTableComponent,
    SearchResultComponent,
    ZernikeDescriptorsListComponent,
    BenchmarkComponent,
    EmdbIdListComponent,
    UploadEmdbIdListFileComponent,
    BenchmarkResultsComponent,
    HomeComponent,
    ReferenceComponent
  ],
  imports: [BrowserModule, HttpClientModule, ReactiveFormsModule, routes],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule {}
