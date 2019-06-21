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
  minutesLeft;
  secondsLeft;
  timeout = false;
  zero = false;
  interval;

  startTimer(descriptor: DescriptorsList) {
    const diffTime =
      Math.abs(
        new Date(descriptor.expirationDate).getTime() - new Date().getTime()
      ) / 1000;
    this.minutesLeft = Math.floor(diffTime / 60);
    this.secondsLeft = Math.floor(diffTime - this.minutesLeft * 60);
    this.interval = setInterval(() => {
      if (this.secondsLeft > 0 || this.minutesLeft > 0) {
        if (this.secondsLeft === 0) {
          this.zero = false;
          this.minutesLeft--;
          this.secondsLeft = 59;
        } else if (this.secondsLeft < 10) {
          this.zero = true;
          this.secondsLeft--;
        } else {
          this.secondsLeft--;
        }
      } else {
        this.timeout = true;
      }
    }, 1000);
  }

  ngOnInit() {
    const emdbIdList = this.route.snapshot.queryParamMap.get("emdbList");
    const contourRepresentation = this.route.snapshot.queryParamMap.get(
      "contour"
    );
    this.descriptorService
      .getDescriptorsList(emdbIdList, parseInt(contourRepresentation))
      .then((data: DescriptorsList) => {
        this.results = data;
        this.showList = true;
        this.startTimer(data);
      });
  }
}
