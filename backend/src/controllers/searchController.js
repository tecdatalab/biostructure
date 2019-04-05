const biomolecule = require("../models/biomolecule");
const sequelize = require("../database").sequelize;

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
      res.status(500);
    } else {
      res.status(200).json(biomolecules);
    }
  } catch (err) {
    res.status(204).json({
      msg: "Biomolecule not found",
      details: err
    });
  }
};

exports.searchResult = async (req, res, next) => {
  emdbid = req.params.emdbID;
  isVolumeFilterOn = req.params.isVolumeFilterOn;
  minRes = req.params.minRes;
  maxRes = req.params.maxRes;
  try {
    let array_results = await getBiomolecules();
    let result = {
      path: "/results/result.hit",
      results: array_results
    };
    res.status(200).json(result);
  } catch (error) {
    res.status(500).json({
      msg: "Biomolecule not found",
      details: err
    });
  }
};

exports.zernike = async (req, res, next) => {
  response = [];
  for (let i = 0; i < 120; i++) {
    response.push(Math.random() * 20);
  }
  res.status(200).json(response);
};

async function getBiomolecules() {
  try {
    let biomolecules = await biomolecule.findAll({
      limit: 8,
      order: [sequelize.random()]
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
