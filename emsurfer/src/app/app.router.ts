import { ModuleWithProviders } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';
import { SearchFormComponent } from './components/search-form/search-form.component';
import { SearchResultComponent } from './components/search-result/search-result.component';

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
    path: 'result/:emdbId',
    component: SearchResultComponent
  },
  {
    path: 'result/emMap',
    component: SearchResultComponent
  }
];

export const routes: ModuleWithProviders = RouterModule.forRoot(router, {onSameUrlNavigation: 'reload'});
