const descriptor = require("../models/descriptorModel");
const fs = require("fs");
var mkdirp = require("mkdirp");
const zipFolder = require("zip-a-folder");
const sequelize = require("../database").sequelize;

exports.getZernikeDescriptors = async (req, res, next) => {
  try {
    const emdbid = req.params.emdbID;
    let descriptorZernike = await descriptor.findOne({
      where: {
        emd_entry_id: parseInt(emdbid)
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
    mkdirp("./public/descriptors/results", function(err) {
      asyncForEach(list, "./public/descriptors/results/", function(response) {
        res.status(200).json(response);
      });
    });
  } catch (err) {
    console.log(err);
    res.status(500).send({
      message: "Backend error"
    });
  }
};

async function asyncForEach(array, path, callback) {
  let zernikeResults = [];
  for (let index = 0; index < array.length; index++) {
    descriptorZernike = await descriptor.findOne({
      where: {
        emd_entry_id: parseInt(array[index])
      }
    });
    if (descriptorZernike) {
      fs.writeFile(
        path + descriptorZernike.emd_entry_id + ".txt",
        descriptorZernike.numbers,
        function(err, result) {
          if (err) console.log("error", err);
        }
      );
      zernikeResults.push(descriptorZernike);
    }
  }
  callback(zernikeResults);
}
