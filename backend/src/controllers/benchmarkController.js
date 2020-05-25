const benchmark_history = require("../models/benchmarkModel");
const biomolecule_emd = require("../models/biomoleculeModelEMD");
const sequelize = require("../database").sequelize;
const fs = require("fs");
var mkdirp = require("mkdirp");
const zipFolder = require("zip-a-folder");

exports.batchQuery = async (req, res) => {
  try {
    const list = req.params.emdblist.split(",");
    let text_info = await getText(req);
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
        mkdirp(benchmarkPath + "/results", function (err) {
           generateFiles(text_info, list, benchmarkPath, newBenchmark.id, function (
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
      })
      .catch(function (err) {
        console.log(err, req.body);
      });
  } catch (err) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

function generateFiles(text_input, IDList, benchmarkPath, idFolder, callback) {
  let final_text = [];
  for (const item of text_input){
    let rank = 1;
    let text = "Rank	EMDB_ID	EUC_D	RESOLUTION \r\n ";
    for (const molecule of item){
      let emb_id = molecule["biomolecule"]["id"];
      let distance = molecule["euc_distance"];
      let resolution = molecule["biomolecule"]["resolution"];
      let temp_txt = rank + "\t" + emb_id + "\t" + distance + "\t" + resolution + "\r\n ";
      text += temp_txt;
      rank += 1;
    }
    final_text.push(text);
  }
  let filesPaths = [];
  let file;
  IDList.forEach(function (ID, index, array) {
    /* Iterate over each ID */
    file = {
      filename: "EMDB-" + ID + ".hit",
      path: benchmarkPath + "/results/EMDB-" + ID + ".hit"
    };
    fs.writeFile(file.path, final_text[index], function (err) {
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
        function (err) {
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

async function getText(req){
  let result_info = []
  try {
    // Query to DB
    let list = req.params.emdblist.split(",");
    for (const ID of list){
      let biomolecules = await sequelize.query(
        "SELECT * FROM top_distance(:emd_id, :type_des, :top)",
        {
          replacements: {
            emd_id: ID,
            type_des: req.params.contour,
            top: req.params.top
          }
        }
      );
      // final result array
      let resultArray = [];
      // Clasification Process
      for (const biomoleculeItem of biomolecules[0]) {
        let bioInfo = await searchByID(biomoleculeItem["emd_id"])
        let distance = biomoleculeItem["distance"].toString();
        resultArray.push({
        biomolecule: bioInfo,
        euc_distance: distance.substring(0, 5)
      })}
      result_info.push(resultArray)
    }
    return result_info;
  } catch (err) {
    return err;
  }
}

async function searchByID(emdbid) {
  return biomolecule_emd
    .findOne({
      where: {
        id: parseInt(emdbid)
      }
    })
    .then(biomolecules => {
      if (!biomolecules) {
        console.log("Biomolecule " + emdbid + " not found.");
        return -1;
      }
      let info = biomolecules["dataValues"];
      return info;
    })
    .catch(err => {
      return -1;
    });
}