const router = require("express").Router();
const uploadFileController = require("../controllers/uploadFileController");

router.post("/EmMap", uploadFileController.uploadFileEmMap2);

module.exports = router;
