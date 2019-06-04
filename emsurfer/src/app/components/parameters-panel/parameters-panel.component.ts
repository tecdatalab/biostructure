import { Component, OnInit } from "@angular/core";
import { ParametersService } from "../../services/parameters.service";
import { Parameters } from "../../models/parameters";
import { UpdateService } from "../../services/update.service";
import { UserService } from "src/app/services/user.service";
import { Update } from "../../models/update";
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
  lastUpdate: Date;

  constructor(
    private userService: UserService,
    private parameterService: ParametersService,
    private updateService: UpdateService,
    private router: Router
  ) {}

  forceUpdate() {
    this.updateService.forceUpdate().then();
  }
  
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
          this.updateService.getLastUpdate().then((response: Update) => {
            this.lastUpdate = response.last_update;
          });
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
