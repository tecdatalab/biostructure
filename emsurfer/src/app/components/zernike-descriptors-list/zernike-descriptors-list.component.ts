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
    const el = document.createElement("textarea");
    el.value = this.values.toString();
    el.value = el.value.replace(/,/g, "\n");
    el.setAttribute("readonly", "");
    document.body.appendChild(el);
    el.select();
    document.execCommand("copy");
    document.body.removeChild(el);
  }
}
