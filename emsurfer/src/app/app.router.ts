import { ModuleWithProviders } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';
import { HomeComponent } from './components/home/home.component';

export const router: Routes = [
  { path: '', redirectTo: '', pathMatch: 'full'}
];

export const routes: ModuleWithProviders = RouterModule.forRoot(router);
