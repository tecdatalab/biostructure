import { Component, OnInit } from "@angular/core";
import { ParametersService } from "../../services/parameters.service";
import { Parameters } from "../../models/parameters";
import { UserService } from "src/app/services/user.service";
import { Router } from "@angular/router";

@Component({
  selector: "app-parameters-panel",
  templateUrl: "./parameters-panel.component.html",
  styleUrls: ["./parameters-panel.component.css"]
})
export class ParametersPanelComponent implements OnInit {
  parameters: Parameters;
  min_volume: number;
  max_volume: number;
  hits_number: number;
  update_rate: number;

  constructor(
    private userService: UserService,
    private router: Router,
    private parameterService: ParametersService
  ) {}

  ngOnInit() {
    if (this.userService.isUserLoggedIn()) {
      this.userService.checkAdminRole().then((data: boolean) => {
        if (data) {
          this.parameterService.getParameters().then((response: Parameters) => {
            this.min_volume = response.volume_filter_min;
            this.max_volume = response.volume_filter_max;
            this.hits_number = response.hits_number;
            this.update_rate = response.update_rate;
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
