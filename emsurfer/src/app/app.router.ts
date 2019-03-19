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
    path:
      'result/:emdb_id/:countour_representation/:volume_filter/:min_res/:max_res',
    component: SearchResultComponent
  },
  {
    path:
      'result/:filename/:countour_level/:countour_representation/:volume_filter/:min_res/:max_res',
    component: SearchResultComponent
  },
  {
    path: 'benchmark',
    component: BenchmarkComponent
  }
];

export const routes: ModuleWithProviders = RouterModule.forRoot(router);
