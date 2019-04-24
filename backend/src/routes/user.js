const router = require("express").Router();
const userController = require("../controllers/userController");

router.post("/auth/token", userController.sendAuthToken);
router.put(
  "/admin/changeUserRole",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  userController.grantAdminRole
);
router.get("/roles", userController.getUsersRoles);
router.get("/users", userController.getUsers);
router.get(
  "/checkAdminRole",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  userController.sendTrueResponse
);

module.exports = router;
