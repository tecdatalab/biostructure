const benchmark = require("../models/benchmarkModel");
const fs = require("fs");
var mkdirp = require("mkdirp");
const zipFolder = require("zip-a-folder");

exports.batchQuery = async (req, res, next) => {
  try {
    await benchmark
      .build({
        date: new Date(),
        fk_user: 0
      })
      .save()
      .then(newBenchmark => {
        const list = req.params.emdblist.split(",");
        const benchmarkPath = "./public/benchmarks/" + newBenchmark.id;
        mkdirp(benchmarkPath + "/results", function(err) {
          generateFiles(list, benchmarkPath, newBenchmark.id, function(
            response
          ) {
            if (Array.isArray(response)) {
              const resJSON = {
                path: benchmarkPath,
                zipFile: "/benchmarks/" + newBenchmark.id + "/results.zip",
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

function generateFiles(IDList, benchmarkPath, idFolder, callback) {
  const text = "Rank	EMDB_ID	EUC_D	RESOLUTION";
  let filesPaths = [];
  let file;
  IDList.forEach(function(ID, index, array) {
    /* Iterate over each ID */
    file = {
      filename: "EMDB-" + ID + ".hit",
      path: benchmarkPath + "/results/EMDB-" + ID + ".hit"
    };
    fs.writeFile(file.path, text, function(err) {
      if (err) {
        callback(err);
      }
    });
    file.path = "/benchmarks/" + idFolder + "/results/EMDB-" + ID + ".hit";
    filesPaths.push(file);
    if (index == array.length - 1) {
      /* Create the zip file */
      zipFolder.zipFolder(
        benchmarkPath + "/results",
        benchmarkPath + "/results.zip",
        function(err) {
          if (err) {
            callback(err);
          } else {
            callback(filesPaths);
          }
        }
      );
    }
  });
}
