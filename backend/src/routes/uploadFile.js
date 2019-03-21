const router = require("express").Router();
const uploadFileController = require("../controllers/uploadFileController");

router.get("/", (req, res, next) => {
    res.send("upload file");
});

router.post(
    '/EmMap',
    uploadFileController.uploadFileEmMap
  );

module.exports = router;
