const router = require("express").Router();
const userController = require("../controllers/userController");

router.post("/auth/token", userController.sendAuthToken);

module.exports = router;
