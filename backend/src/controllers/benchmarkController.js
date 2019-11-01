const benchmark_history = require("../models/benchmarkModel");
const fs = require("fs");
var mkdirp = require("mkdirp");
const zipFolder = require("zip-a-folder");

exports.batchQuery = async (req, res, next) => {
  try {
    const list = req.params.emdblist.split(",");
    await benchmark_history
      .build({
        date_time: new Date(),
        ip: req.connection.remoteAddress,
        user_id: null,
        representation_id: req.params.contour,
        volume_filter_id: req.params.filter == "true" ? 1 : 0,
        top_results: req.params.top,
        emd_list: list
      })
      .save()
      .then(newBenchmark => {
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
              console.log("I/O Error");
              res.status(500).send({
                message: "I/O Server error"
              });
            }
          });
        });
      });
  } catch (err) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

function generateFiles(IDList, benchmarkPath, idFolder, callback) {
  const text =
    "Rank	EMDB_ID	EUC_D	RESOLUTION \r\n 1	6409	12.851	22 \r\n 2	4804	12.945	14 \r\n 3	6478	15.250	35 \r\n 4	9618	19.548	10.6 \r\n";
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
