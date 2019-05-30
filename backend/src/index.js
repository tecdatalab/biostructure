const express = require("express");
const app = express();
const cors = require("cors");
const bodyParser = require("body-parser");
const expressValidator = require("express-validator");
const indexRoutes = require("./routes/index");
const searchRoutes = require("./routes/search");
const uploadFileRoutes = require("./routes/uploadFile");
const benchmarkRoutes = require("./routes/benchmark");
const contourRepresentationRoutes = require("./routes/contourRepresentation");
const checkerRoutes = require("./routes/checker");
const userRoutes = require("./routes/user");
const descriptorsRoutes = require("./routes/descriptor");
const searchHistoryRoutes = require("./routes/searchHistory");
const parameters = require("./routes/parameters");
//settings
app.set("port", process.env.PORT || 3000);

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
app.use("/contour", contourRepresentationRoutes);
app.use("/checker", checkerRoutes);
app.use("/user", userRoutes);
app.use("/descriptor", descriptorsRoutes);
app.use("/history", searchHistoryRoutes);
app.use("/parameters", parameters);

//public routes
app.use("/img", express.static("public/img"));
app.use("/results", express.static("public/results"));
app.use("/benchmarks", express.static("public/benchmarks"));
app.use("/descriptors", express.static("public/descriptors"));

app.listen(app.get("port"), () => {
  console.log("Server on port ", app.get("port"));
});
