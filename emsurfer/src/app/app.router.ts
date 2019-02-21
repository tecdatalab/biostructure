import { ModuleWithProviders } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

export const router: Routes = [
  { path: '', redirectTo: '', pathMatch: 'full'}
];

export const routes: ModuleWithProviders = RouterModule.forRoot(router);
