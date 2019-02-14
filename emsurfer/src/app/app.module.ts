import { BrowserModule } from '@angular/platform-browser';
import { RouterModule, Routes } from '@angular/router';
import { NgModule } from '@angular/core';
import { ReactiveFormsModule } from '@angular/forms';

import { routes } from './app.router';

import { AppComponent } from './app.component';
import { HeaderComponent} from './components/header/header.component'
import { ContourShapeInputComponent } from './components/contour-shape-input/contour-shape-input.component';
import { VolumeFilterInputComponent } from './components/volume-filter-input/volume-filter-input.component';

const appRoutes: Routes = [
  {
    path: '', 
    redirectTo: 'search', 
    pathMatch: 'full'
  },
  {
    path: 'search', 
    component: VolumeFilterInputComponent
  }
]

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    ContourShapeInputComponent,
    VolumeFilterInputComponent
  ],
  imports: [
    BrowserModule,
    ReactiveFormsModule,
    RouterModule.forRoot(
      appRoutes
    )
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
