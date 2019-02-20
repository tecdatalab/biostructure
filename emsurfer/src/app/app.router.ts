import { ModuleWithProviders } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';
import { ContourShapeInputComponent } from './components/contour-shape-input/contour-shape-input.component';

export const router: Routes = [
  {
    path: '',
    redirectTo: 'search',
    pathMatch: 'full'
  },
  {
    path: 'search',
    component: ContourShapeInputComponent
  }
];

export const routes: ModuleWithProviders = RouterModule.forRoot(router);
