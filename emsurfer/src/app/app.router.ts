import { ModuleWithProviders } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';
import { SearchFormComponent } from './components/search-form/search-form.component';
import { SearchResultComponent } from './components/search-result/search-result.component';
import { BenchmarkComponent } from './components/benchmark/benchmark.component';

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
  },
  {
    path: 'benchmark',
    component: BenchmarkComponent
  }
];

export const routes: ModuleWithProviders = RouterModule.forRoot(router);
