import { Component, OnInit, Input } from "@angular/core";
import { Descriptor } from "src/app/models/descriptor";

@Component({
  selector: "app-zernike-result",
  templateUrl: "./zernike-result.component.html",
  styleUrls: ["./zernike-result.component.css"]
})
export class ZernikeResultComponent implements OnInit {
  @Input() zernike: Descriptor;
  validZernike = true;
  constructor() {}

  ngOnInit() {
    if (!this.zernike.numbers) {
      this.validZernike = false;
      console.log(this.zernike.numbers);
    }
  }
}
