import { BrowserModule } from '@angular/platform-browser';
import { NgModule } from '@angular/core';
import { ReactiveFormsModule } from '@angular/forms';

import { routes } from './app.router';

import { AppComponent } from './app.component';
import { HeaderComponent } from './components/header/header.component';
import { ContourShapeInputComponent } from './components/contour-shape-input/contour-shape-input.component';
import { VolumeFilterInputComponent } from './components/volume-filter-input/volume-filter-input.component';
import { SearchFormComponent } from './components/search-form/search-form.component';
import { HomeComponent } from './components/home/home.component';
import { EmdbIdInputComponent } from './components/emdb-id-input/emdb-id-input.component';
import { UploadEmMapComponent } from './components/upload-em-map/upload-em-map.component';
import { QueryMethodComponent } from './components/query-method/query-method.component';
import { ResolutionFilterComponent } from './components/resolution-filter/resolution-filter.component';
import { BiomoleculeComponent } from './components/biomolecule/biomolecule.component';
import { BiomoleculesTableComponent } from './components/biomolecules-table/biomolecules-table.component';
import { SearchResultComponent } from './components/search-result/search-result.component';
import { ZernikeDescriptorsListComponent } from './components/zernike-descriptors-list/zernike-descriptors-list.component';

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    ContourShapeInputComponent,
    VolumeFilterInputComponent,
    SearchFormComponent,
    HomeComponent,
    EmdbIdInputComponent,
    UploadEmMapComponent,
    QueryMethodComponent,
    ResolutionFilterComponent,
    BiomoleculeComponent,
    BiomoleculesTableComponent,
    SearchResultComponent,
    ZernikeDescriptorsListComponent
  ],
  imports: [BrowserModule, ReactiveFormsModule, routes],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule {}
