const biomolecule = require("../models/biomoleculeModel");
const biomolecule_pdb = require("../models/biomoleculeModelPDB");
const search_history = require("../models/searchHistoryModel");
const Op = require("../database").Op;

exports.searchByID = async (req, res) => {
  try {
    const table = req.params.app;
    if (parseInt(table)==1){
      const emdbid = req.params.ID;
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
    } else {
      const pdbid = req.params.ID;
      let biomolecules = await biomolecule_pdb.findOne({
        where: {
          pdb_id: pdbid
        }
      });
      if (!biomolecules) {
        console.log("Biomolecule " + pdbid + " not found.");
        res.status(400).send({
          message: "Biomolecule " + pdbid + " not found."
        });
      } else {
        res.status(200).json(biomolecules);
      }
    }
  } catch (err) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

exports.saveSearch = async (req, res, next) => {
  try {
    await search_history
      .build({
        date_time: new Date(),
        ip: req.connection.remoteAddress,
        user_id: req.authenticatedUser.id,
        emd_entry_id: req.params.emdbID == null ? 1 : req.params.emdbID,
        name_file: req.params.filename,
        contour_level: req.params.contourLevel,
        representation_id: req.params.contourRepresentation,
        volume_filter_id: req.params.isVolumeFilterOn == "true" ? 1 : 0,
        resolution_filter_min: checkMinResolutionFilter(req.params.minRes),
        resolution_filter_max: checkMaxResolutionFilter(req.params.maxRes)
      })
      .save()
      .then(next());
  } catch (err) {
    res.status(500).send({
      message: "Internal server error"
    });
  }
};

exports.searchResult = async (req, res) => {
  try {
    minRes = checkMinResolutionFilter(req.params.minRes);
    maxRes = checkMaxResolutionFilter(req.params.maxRes);
    let query_results = await getBiomolecules(minRes, maxRes);
    let result = {
      path: "/results/result.hit",
      results: query_results
    };
    res.status(200).json(result);
  } catch (error) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

exports.searchResultMap = async (req, res) => {
  try {
    minRes = checkMinResolutionFilter(req.params.minRes);
    maxRes = checkMaxResolutionFilter(req.params.maxRes);
    let query_results = await getBiomolecules(minRes, maxRes);
    let result = {
      path: "/results/result.hit",
      results: query_results
    };
    res.status(200).json(result);
  } catch (error) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

exports.getCathDetail = async (req, res) => {
  try {
    const pdb_id = req.params.ID;
    let cath = await cathInfo.findOne({
      where: {
        atomic_structure_id: pdb_id
      }
    });
    if (!cath) {
      console.log("Cath information " + pdb_id + " not found.");
      res.status(400).send({
        message: "Cath information " + pdb_id + " not found."
      });
    } else {
      res.status(200).json(cath);
    }
  }catch (error) {
    res.status(500).send({
      message: "Backend error"
    });
  }
}

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
    Search into the database to set the default value of min resolution.
    */
    return 2;
  } else {
    return parseFloat(minRes);
  }
}

function checkMaxResolutionFilter(maxRes) {
  if (isNaN(maxRes)) {
    /*
    Search into the database to set the default value of max resolution.
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
        euc_distance: 5
      });
    });
    return resultArray;
  } catch (err) {
    return err;
  }
}
