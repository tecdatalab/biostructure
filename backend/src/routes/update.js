const router = require("express").Router();
const updateController = require("../controllers/updateController");

router.get("/", updateController.getLastUpdate);

module.exports = router;
