const biomolecule = require("../models/biomoleculeModel");
const sequelize = require("../database").sequelize;
const Op = require("../database").Op;

exports.searchByID = async (req, res, next) => {
  emdbid = req.params.emdbID;
  try {
    let biomolecules = await biomolecule.findOne({
      where: {
        id: parseInt(emdbid)
      }
    });
    if (!biomolecules) {
      console.log("Biomolecule " + emdbid + " not found.");
      res.status(400).send({
        message: "Biomolecule " + emdbid + " not found."
      });
    } else {
      res.status(200).json(biomolecules);
    }
  } catch (err) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

exports.searchResult = async (req, res, next) => {
  emdbid = req.params.emdbID;
  contourRepresentation = req.params.contourRepresentation;
  isVolumeFilterOn = req.params.isVolumeFilterOn;
  minRes = req.params.minRes;
  maxRes = req.params.maxRes;
  try {
    let array_results = await getBiomolecules(minRes, maxRes);
    let result = {
      path: "/results/result.hit",
      results: array_results
    };
    res.status(200).json(result);
  } catch (error) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

exports.searchResultMap = async (req, res, next) => {
  emdbid = req.params.emdbID;
  contourRepresentation = req.params.contourRepresentation;
  isVolumeFilterOn = req.params.isVolumeFilterOn;
  minRes = req.params.minRes;
  maxRes = req.params.maxRes;
  try {
    let array_results = await getBiomolecules(minRes, maxRes);
    let result = {
      path: "/results/result.hit",
      results: array_results
    };
    res.status(200).json(result);
  } catch (error) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

exports.zernike = async (req, res, next) => {
  response = [];
  for (let i = 0; i < 120; i++) {
    response.push((Math.random() * 20).toFixed(3));
  }
  res.status(200).json(response);
};

exports.zernikeMap = async (req, res, next) => {
  response = [];
  for (let i = 0; i < 120; i++) {
    response.push(Math.random() * 20);
  }
  res.status(200).json(response);
};

async function getBiomolecules(minRes, maxRes) {
  console.log("Min res: " + minRes);
  console.log("Max res: " + maxRes);
  try {
    let biomolecules = await biomolecule.findAll({
      where: {
        volume: {
          [Op.gte]: minRes,
          [Op.lte]: maxRes
        }
      }
    });
    let resultArray = [];
    biomolecules.forEach(biomoleculeItem => {
      resultArray.push({
        biomolecule: biomoleculeItem.get({
          plain: true
        }),
        euc_distance: 5,
        ratio_of_volume: 4,
        resolution: 3
      });
    });
    return resultArray;
  } catch (err) {
    return err;
  }
}
