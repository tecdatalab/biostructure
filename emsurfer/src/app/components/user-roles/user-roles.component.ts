import { Component, OnInit } from "@angular/core";
import { User } from "src/app/models/user";

@Component({
  selector: "app-user-roles",
  templateUrl: "./user-roles.component.html",
  styleUrls: ["./user-roles.component.css"]
})
export class UserRolesComponent implements OnInit {
  constructor() {}
  users = [];
  filter;
  roles = [2, 0, 1];
  verify(user) {
    if (this.filter) {
      if (user.id == this.filter) {
        return true;
      } else {
        return false;
      }
    }
    return true;
  }

  createValues() {
    for (let i = 0; i < 10000; i++) {
      this.users.push({
        id: i,
        name: "Name" + i,
        email: i + "@asdf.com",
        role: i % 2
      });
    }
  }

  ngOnInit() {
    this.createValues();
  }
}
