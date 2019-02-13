import { BrowserModule } from '@angular/platform-browser';
import { RouterModule, Routes } from '@angular/router';
import { NgModule } from '@angular/core';
import { FormsModule, ReactiveFormsModule } from '@angular/forms';

import { AppComponent } from './app.component';
import { HeaderComponent} from './components/header/header.component'
import { HomeComponent } from './components/home/home.component'
import { ContourShapeInputComponent } from './components/contour-shape-input/contour-shape-input.component';
import { UploadEmMapComponent } from './components/upload-em-map/upload-em-map.component';
import { SearchFormComponent } from './components/search-form/search-form.component';

const appRoutes: Routes = [
  {
    path: '', 
    redirectTo: 'search', 
    pathMatch: 'full'
  },
  {
    path: 'home', 
    component: ContourShapeInputComponent
  },
  {
    path: 'search', 
    component: SearchFormComponent
  },
  {
    path: 'query/:id/:volume/:min/:max',
    component: HomeComponent
  }
]

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    HomeComponent,
    ContourShapeInputComponent,
    UploadEmMapComponent,
    SearchFormComponent
  ],
  imports: [
    BrowserModule,
    ReactiveFormsModule,
    FormsModule,
    RouterModule.forRoot(
      appRoutes
    )
  ],
  exports: [
    RouterModule,
    HeaderComponent,
    HomeComponent,
    ContourShapeInputComponent,
    UploadEmMapComponent,
    SearchFormComponent
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
