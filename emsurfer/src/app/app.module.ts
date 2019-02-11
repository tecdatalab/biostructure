import { BrowserModule } from '@angular/platform-browser';
import { RouterModule, Routes } from '@angular/router';
import { NgModule } from '@angular/core';
import { ReactiveFormsModule } from '@angular/forms';

import { HeaderComponent} from './components/header/header.component'
import { AppComponent } from './app.component';
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
  }
]

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    ContourShapeInputComponent,
    UploadEmMapComponent,
    SearchFormComponent
  ],
  imports: [
    BrowserModule,
    ReactiveFormsModule,
    RouterModule.forRoot(
      appRoutes,
      {useHash: true}
    )
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
