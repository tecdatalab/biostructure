import { ModuleWithProviders } from "@angular/core";
import { Routes, RouterModule } from "@angular/router";
import { SearchFormComponent } from "./components/search-form/search-form.component";
import { SearchResultComponent } from "./components/search-result/search-result.component";
import { BenchmarkComponent } from "./components/benchmark/benchmark.component";
import { BenchmarkResultsComponent } from "./components/benchmark-results/benchmark-results.component";
import { HomeComponent } from "./components/home/home.component";

export const router: Routes = [
  {
    path: "",
    redirectTo: "home",
    pathMatch: "full"
  },
  {
    path: "search",
    component: SearchFormComponent
  },
  {
    path: "result/:emdbId",
    component: SearchResultComponent
  },
  {
    path: "result/emMap",
    component: SearchResultComponent
  },
  {
    path: "benchmark",
    component: BenchmarkComponent
  },
  {
    path: "benchmark/results",
    component: BenchmarkResultsComponent
  },
  {
    path: "home",
    component: HomeComponent
  }
];

export const routes: ModuleWithProviders = RouterModule.forRoot(router);
