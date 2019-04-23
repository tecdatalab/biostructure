import { Component, OnInit } from "@angular/core";
import { User } from "src/app/models/user";
import { UserService } from "src/app/services/user.service";

@Component({
  selector: "app-user-roles",
  templateUrl: "./user-roles.component.html",
  styleUrls: ["./user-roles.component.css"]
})
export class UserRolesComponent implements OnInit {
  constructor() {}
  users = [];
  checkedOption;
  currentPage = -1;
  filter;
  roles = [2, 0, 1];

  verify(user) {
    if (this.filter) {
      if (user[this.checkedOption].includes(this.filter)) {
        this.currentPage = -1;
        return true;
      } else {
        return false;
      }
    }
    return true;
  }

  createValues() {
    for (let i = 0; i < 150; i++) {
      this.users.push({
        id: i,
        name: "Name" + i,
        email: i + "@asdf.com",
        role: i % 2
      });
    }
  }

  nextPage() {
    if (this.currentPage + 100 < this.users.length) {
      this.currentPage += 100;
    }
  }

  previousPage() {
    if (this.currentPage > 0) {
      this.currentPage -= 100;
    }
  }

  ngOnInit() {
    this.createValues();
  }
}
