const descriptor = require("../models/descriptorModel");
const fs = require("fs");
var mkdirp = require("mkdirp");
const zipFolder = require("zip-a-folder");
const createIndex = require("../utilities/indexGenerator");

exports.getZernikeDescriptors = async (req, res, next) => {
  try {
    const emdbid = req.params.emdbID;
    const type_descriptor = req.params.contourRepresentation;
    let descriptorZernike = await descriptor.findOne({
      where: {
        emd_entry_id: parseInt(emdbid),
        type_descriptor_id: parseInt(type_descriptor)
      }
    });
    if (!descriptorZernike) {
      console.log("Biomolecule descriptors " + emdbid + " not found.");
      res.status(400).send({
        message: "Biomolecule descriptors " + emdbid + " not found."
      });
    } else {
      res.status(200).json(descriptorZernike);
    }
  } catch (err) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

exports.getZernikeList = async (req, res, next) => {
  try {
    const list = req.params.emdblist.split(",");
    const type_descriptor = req.params.contourRepresentation;
    const pathFiles = "./public/descriptors/" + createIndex.createIndex();
    mkdirp(pathFiles + "/results", function(err) {
      asyncForEach(list, pathFiles, type_descriptor, function(response) {
        if (response instanceof Error) {
          console.log("I/O Error");
          res.status(500).send({
            message: "I/O Server error"
          });
        } else {
          res.status(200).json(response);
        }
      });
    });
  } catch (err) {
    console.log(err);
    res.status(500).send({
      message: "Backend error"
    });
  }
};

async function asyncForEach(array, path, contour, callback) {
  let zernikeResults = [];
  //generate each zernike descriptor file
  for (let index = 0; index < array.length; index++) {
    descriptorZernike = await descriptor.findOne({
      where: {
        emd_entry_id: parseInt(array[index]),
        type_descriptor_id: contour
      }
    });
    if (descriptorZernike) {
      fs.writeFile(
        path + "/results/" + descriptorZernike.emd_entry_id + ".txt",
        descriptorZernike.numbers,
        function(err, result) {
          if (err) console.log("error", err);
        }
      );
      zernikeResults.push(descriptorZernike);
    } else {
      zernikeResults.push({
        emd_entry_id: array[index],
        type_descriptor_id: 0
      });
    }
  }
  //generate the compressed file
  zipFolder.zipFolder(path + "/results", path + "/results.zip", function(err) {
    if (err) {
      callback(err);
    } else {
      var time = new Date();
      time.setMinutes(time.getMinutes() + 10);
      callback({
        path: path.substring(8) + "/results.zip",
        expirationDate: time,
        descriptors: zernikeResults
      });
    }
  });
}
