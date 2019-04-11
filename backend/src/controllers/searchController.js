const biomolecule = require("../models/biomoleculeModel");
const search_history = require("../models/searchHistoryModel");
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
    minRes = checkMinResolutionFilter(minRes);
    maxRes = checkMaxResolutionFilter(maxRes);
    let array_results = await getBiomolecules(minRes, maxRes);
    await search_history
      .build({
        date_time: new Date(),
        ip: req.connection.remoteAddress,
        emd_entry_id: emdbid,
        name_file: "-",
        counter_level: 0.0,
        representation_id: contourRepresentation,
        volume_filter_id: 0,
        resolution_filter_min: minRes,
        resolution_filter_max: maxRes
      })
      .save()
      .then(newSearch => {
        let result = {
          path: "/results/result.hit",
          results: array_results
        };
        res.status(200).json(result);
      });
  } catch (error) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

exports.searchResultMap = async (req, res, next) => {
  filename = req.params.filename;
  contourRepresentation = req.params.contourRepresentation;
  contourLevel = req.params.contourLevel;
  isVolumeFilterOn = req.params.isVolumeFilterOn;
  minRes = req.params.minRes;
  maxRes = req.params.maxRes;
  try {
    minRes = checkMinResolutionFilter(minRes);
    maxRes = checkMaxResolutionFilter(maxRes);
    let array_results = await getBiomolecules(minRes, maxRes);
    await search_history
      .build({
        date_time: new Date(),
        ip: req.connection.remoteAddress,
        emd_entry_id: 1,
        name_file: filename,
        counter_level: contourLevel,
        representation_id: contourRepresentation,
        volume_filter_id: 0,
        resolution_filter_min: minRes,
        resolution_filter_max: maxRes
      })
      .save()
      .then(newSearch => {
        let result = {
          path: "/results/result.hit",
          results: array_results
        };
        res.status(200).json(result);
      })
      .catch(err => {
        console.log(err);
      });
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

function checkMinResolutionFilter(minRes) {
  if (isNaN(minRes)) {
    /*
    Search in the database to set the default value of min resolution.
    */
    return 2;
  } else {
    return parseFloat(minRes);
  }
}

function checkMaxResolutionFilter(maxRes) {
  if (isNaN(maxRes)) {
    /*
    Search in the database to set the default value of max resolution.
    */
    return 10;
  } else {
    return parseFloat(maxRes);
  }
}

async function getBiomolecules(minRes, maxRes) {
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
