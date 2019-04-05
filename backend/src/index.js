const express = require("express");
const app = express();
const cors = require("cors");
const bodyParser = require("body-parser");
const expressValidator = require("express-validator");
const indexRoutes = require("./routes/index");
const searchRoutes = require("./routes/search");
const uploadFileRoutes = require("./routes/uploadFile");
const benchmarkRoutes = require("./routes/benchmark");

//settings
app.set("port", process.env.PORT || 3001);

//middlewares
app.use(bodyParser.urlencoded({ limit: "50mb", extended: true }));
app.use(bodyParser.json({ limit: "50mb", extended: true }));
app.use(express.json());
app.use(cors());
app.use(expressValidator());

//routes
app.use(indexRoutes);
app.use("/search", searchRoutes);
app.use("/upload", uploadFileRoutes);
app.use("/benchmark", benchmarkRoutes);

//public routes
app.use("/img", express.static("public/img"));
app.use("/results", express.static("public/results"));
app.use("/benchmarks", express.static("public/benchmarks"));

app.listen(app.get("port"), () => {
  console.log("Server on port ", app.get("port"));
});
