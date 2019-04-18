import { Component, OnInit } from "@angular/core";

@Component({
  selector: "app-zernike-descriptors-module",
  templateUrl: "./zernike-descriptors-module.component.html"
})
export class ZernikeDescriptorsModuleComponent implements OnInit {
  constructor() {}

  values;
  descriptors;
  ngOnInit() {
    this.values = [8, 9, 7, 6];
    this.descriptors = [1, 2, 3, 4];
  }
}
