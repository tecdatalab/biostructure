import { Component, OnInit } from "@angular/core";
import { ActivatedRoute } from "@angular/router";
import { DescriptorService } from "src/app/services/descriptor.service";
import { DescriptorsList } from "src/app/models/descriptorsList";

@Component({
  selector: "app-zernike-descriptors-module",
  templateUrl: "./zernike-descriptors-module.component.html",
  providers: [DescriptorService]
})
export class ZernikeDescriptorsModuleComponent implements OnInit {
  constructor(
    private descriptorService: DescriptorService,
    private route: ActivatedRoute
  ) {}
  results: DescriptorsList;
  showList = false;
  ngOnInit() {
    const emdbList = this.route.snapshot.paramMap.get("emdbList");
    this.descriptorService
      .getDescriptorsList(emdbList)
      .then((data: DescriptorsList) => {
        this.results = data;
        this.showList = true;
      });
  }
}
