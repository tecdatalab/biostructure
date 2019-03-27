import { Component, Input } from "@angular/core";

@Component({
  selector: "app-zernike-descriptors-list",
  templateUrl: "./zernike-descriptors-list.component.html",
  styleUrls: ["zernike-descriptors-list.component.css"]
})
export class ZernikeDescriptorsListComponent {
  @Input() descriptors: number[];
  @Input() values: number[];
  constructor() {}

  cpyToClipboard() {
    const el = document.createElement("input");
    el.value = this.values.toString();
    el.value.replace(",", "\n");
    el.setAttribute("readonly", "");
    document.body.appendChild(el);
    el.select();
    document.execCommand("copy");
    document.body.removeChild(el);
  }
}
