import { ModuleWithProviders } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';
import { SearchFormComponent } from './components/search-form/search-form.component';
import { HomeComponent } from './components/home/home.component';

export const router: Routes = [
  {
    path: '',
    redirectTo: 'search',
    pathMatch: 'full'
  },
  {
    path: 'search',
    component: SearchFormComponent
  },
  {
    path: 'query/:id/:volume/:min/:max',
    component: HomeComponent
  }
];

export const routes: ModuleWithProviders = RouterModule.forRoot(router);
