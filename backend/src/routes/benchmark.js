const router = require("express").Router();
const benchmarkController = require("../controllers/benchmarkController");

router.post("/query", benchmarkController.batchQuery);

module.exports = router;
