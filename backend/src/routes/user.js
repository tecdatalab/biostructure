const router = require("express").Router();
const userController = require("../controllers/userController");

router.post("/auth/token", userController.sendAuthToken);
router.put(
  "/admin/changeUserRole",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  userController.changeUserRole
);
router.get("/roles", userController.getUsersRoles);
router.get("/users", userController.getUsers);
router.post(
  "/checkAdminRole",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  userController.sendTrueResponse
);

module.exports = router;
