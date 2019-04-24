const router = require("express").Router();
const userController = require("../controllers/userController");

router.post("/auth/token", userController.sendAuthToken);
router.put(
  "/admin/grantAdminRole",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  userController.grantAdminRole
);
router.get("/roles", userController.getUsersRoles);
router.get("/users", userController.getUsers);

module.exports = router;
