const router = require("express").Router();
const userController = require("../controllers/userController");

router.post("/auth/token", userController.sendAuthToken);
router.put(
  "/admin/changeUserRole",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  userController.changeUserRole
);
router.get(
  "/roles",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  userController.getUsersRoles
);
router.get(
  "/users",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  userController.getUsers
);
router.get(
  "/checkAdminRole",
  userController.verifyUserToken,
  userController.isUserAdmin
);

module.exports = router;
