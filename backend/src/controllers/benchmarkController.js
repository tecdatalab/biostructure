const benchmark = require("../models/benchmark");
const sequelize = require("../database").sequelize;
const fs = require("fs");
var archiver = require("archiver");
var mkdirp = require("mkdirp");

exports.batchQuery = async (req, res, next) => {
  try {
    await benchmark
      .build({
        date: new Date(),
        fk_user: 0
      })
      .save()
      .then(newBenchmark => {
        console.log("BODY " + req.body);
        const benchmarkPath = "./public/benchmarks/" + newBenchmark.id;
        mkdirp(benchmarkPath, function(err) {
          generateFiles(req.body, benchmarkPath, function(response) {
            if (Array.isArray(response)) {
              const resJSON = {
                path: benchmarkPath,
                zipFile: "",
                results: response
              };
              res.status(200).json(resJSON);
            } else {
              console.log("Error");
              res.status(204).json({
                msg: "ERROR: I/O error.",
                details: response
              });
            }
          });
        });
      });
  } catch (err) {
    res.status(204).json({
      msg: "ERROR",
      details: err
    });
  }
};

function generateFiles(IDList, benchmarkPath, callback) {
  const text = "Rank	EMDB_ID	EUC_D	RESOLUTION";
  let filesPaths = [];
  let file;
  IDList.forEach(function(ID, index, array) {
    file = {
      filename: "EMDB-" + ID + ".hit",
      path: benchmarkPath + "/EMDB-" + ID + ".hit"
    };
    fs.writeFile(file.path, text, function(err) {
      if (err) {
        callback(err);
      }
    });
    file.path = file.path.substring(1);
    filesPaths.push(file);
    if (index == array.length - 1) {
      callback(filesPaths);
    }
  });
}
