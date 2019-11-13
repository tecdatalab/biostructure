const router = require("express").Router();
const benchmarkController = require("../controllers/benchmarkController");

router.get(
  "/query/:emdblist/:contour/:filter/:top",
  benchmarkController.batchQuery
);

module.exports = router;
