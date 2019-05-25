const router = require("express").Router();
const updateController = require("../controllers/updateController");
const userController = require("../controllers/userController");

router.get(
  "/",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  updateController.getLastUpdate
);

router.get(
  "/forceUpdater",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  updateController.forceUpdate
);

module.exports = router;
