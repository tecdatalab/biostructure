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

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    ContourShapeInputComponent,
    VolumeFilterInputComponent,
    SearchFormComponent,
    HomeComponent
  ],
  imports: [BrowserModule, ReactiveFormsModule, routes],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule {}
