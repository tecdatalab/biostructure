const router = require("express").Router();
const parametersController = require("../controllers/parametersController");
const userController = require("../controllers/userController");

router.get(
  "/",
  userController.verifyUserToken,
  userController.verifyAdminUser,
  parametersController.getParameters
);

module.exports = router;
