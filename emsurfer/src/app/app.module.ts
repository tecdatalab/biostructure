import { BrowserModule } from '@angular/platform-browser';
import { NgModule } from '@angular/core';

import { HeaderComponent} from './components/header/header.component'

import { AppComponent } from './app.component';
import { ContourShapeInputComponent } from './components/contour-shape-input/contour-shape-input.component';
import { UploadEmMapComponent } from './components/upload-em-map/upload-em-map.component';
import { SearchFormComponent } from './components/search-form/search-form.component';

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    ContourShapeInputComponent,
    UploadEmMapComponent,
    SearchFormComponent
  ],
  imports: [
    BrowserModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
