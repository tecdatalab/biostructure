const router = require("express").Router();
const parametersController = require("../controllers/parametersController");
const userController = require("../controllers/userController");

router.get(
  "/",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  parametersController.getParameters
);

router.get(
  "/set/:volumeMin/:volumeMax/:hits/:update",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  parametersController.setParameters
);

module.exports = router;
