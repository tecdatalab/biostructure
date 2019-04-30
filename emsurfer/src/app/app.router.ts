import { ModuleWithProviders } from "@angular/core";
import { Routes, RouterModule } from "@angular/router";
import { SearchFormComponent } from "./components/search-form/search-form.component";
import { SearchResultComponent } from "./components/search-result/search-result.component";
import { BenchmarkComponent } from "./components/benchmark/benchmark.component";
import { BenchmarkResultsComponent } from "./components/benchmark-results/benchmark-results.component";
import { HomeComponent } from "./components/home/home.component";
import { ReferenceComponent } from "./components/reference/reference.component";
import { ContactComponent } from "./components/contact/contact.component";

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
  },
  {
    path: "reference",
    component: ReferenceComponent
  },
  {
    path: "contact",
    component: ContactComponent
  }
];

export const routes: ModuleWithProviders = RouterModule.forRoot(router);
