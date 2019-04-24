import { Component, OnInit } from "@angular/core";
import { User } from "src/app/models/user";
import { UserRole } from "src/app/models/userRole";
import { UserService } from "src/app/services/user.service";
import { Router } from "@angular/router";

@Component({
  selector: "app-user-roles",
  templateUrl: "./user-roles.component.html"
})
export class UserRolesComponent implements OnInit {
  constructor(private userService: UserService, private router: Router) {}
  users: User[];
  checkedOption = "name";
  selectedRoleFilter: number;
  currentPage = -1;
  value;
  roles: UserRole[];

  filterFunction(collection) {
    return collection.filter(user => {
      if (this.value) {
        if (
          user[this.checkedOption].includes(this.value) &&
          (user["role"] == this.selectedRoleFilter ||
            this.selectedRoleFilter == -1)
        ) {
          return true;
        }
        return false;
      } else if (
        user["role"] != this.selectedRoleFilter &&
        this.selectedRoleFilter != -1
      ) {
        return false;
      }
      return true;
    });
  }

  onChangeRole(value, user) {
    this.userService.changeUserRole(
      user.id,
      parseInt(this.roles[parseInt(value)].id)
    );
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
    if (this.userService.isUserLoggedIn()) {
      this.userService.checkAdminRole().then((data: boolean) => {
        if (data) {
          console.log("DATA " + data);
          this.selectedRoleFilter = -1;
          this.userService.getUsers().then((data: User[]) => {
            this.users = data;
          });
          this.userService.getUserRoles().then((data: UserRole[]) => {
            this.roles = data;
          });
        }
      });
    } else {
      this.router.navigate(["/home"]);
    }
  }
}
