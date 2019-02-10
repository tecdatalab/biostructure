import { BrowserModule } from '@angular/platform-browser';
import { NgModule } from '@angular/core';

import { HeaderComponent} from './header/header.component'

import { AppComponent } from './app.component';
import { ContourShapeInputComponent } from './contour-shape-input/contour-shape-input.component';
@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    ContourShapeInputComponent
  ],
  imports: [
    BrowserModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
