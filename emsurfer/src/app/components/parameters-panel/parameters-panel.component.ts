import { Component, OnInit } from "@angular/core";
import { ParametersService } from "../../services/parameters.service";
import { Parameters } from "../../models/parameters";
import { UserService } from "src/app/services/user.service";
import { Router } from "@angular/router";

@Component({
  selector: "app-parameters-panel",
  templateUrl: "./parameters-panel.component.html"
})
export class ParametersPanelComponent implements OnInit {
  parameters: Parameters;
  minVolume: number;
  maxVolume: number;
  hitsNumber: number;
  updateRate: number;

  constructor(
    private userService: UserService,
    private router: Router,
    private parameterService: ParametersService
  ) {}

  saveChanges() {
    this.parameterService
      .setParameters(
        this.minVolume,
        this.maxVolume,
        this.hitsNumber,
        this.updateRate
      )
      .then(response => {
        alert("Changes saved successfully");
      });
  }

  ngOnInit() {
    if (this.userService.isUserLoggedIn()) {
      this.userService.checkAdminRole().then((data: boolean) => {
        if (data) {
          this.parameterService.getParameters().then((response: Parameters) => {
            this.minVolume = response.volume_filter_min;
            this.maxVolume = response.volume_filter_max;
            this.hitsNumber = response.hits_number;
            this.updateRate = response.update_rate;
          });
        } else {
          this.router.navigate(["/home"]);
        }
      });
    } else {
      this.router.navigate(["/home"]);
    }
  }
}
