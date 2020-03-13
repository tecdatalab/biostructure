const biomolecule_emd = require("../models/biomoleculeModelEMD");
const biomolecule_pdb = require("../models/biomoleculeModelPDB");
const cathInfo_pdb = require("../models/cathModel");
const sequelize = require("../database").sequelize;
const search_history = require("../models/searchHistoryModel");
const Op = require("../database").Op;

//#region General Functions

exports.searchByID = async (req, res) => {
  try {
    const table = req.params.app;
    let id = req.params.ID;
    let biomolecule = null;
    if (table == "emdb") {
      biomolecule = await biomolecule_emd.findOne({
        where: {
          id: parseInt(id)
        }
      });
    } else {
      biomolecule = await biomolecule_pdb.findOne({
        where: {
          id_code: id
        }
      });
    }
    giveResponse(id, biomolecule, res);
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

//#endregion

//#region Functions for EMDB Queries

exports.searchResult = async (req, res) => {
  try {
    emd_id = req.params.emdbID;
    type_descriptor = req.params.typeDescriptor;
    top_can = req.params.topCan;
    is_volume_filter_on = req.params.isVolumeFilterOn;
    min_res = req.params.minRes;
    max_res = req.params.maxRes;

    let query_results = await getBiomolecules(
      emd_id,
      type_descriptor,
      top_can,
      is_volume_filter_on,
      min_res,
      max_res
    );

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

//#endregion

//#region Functions for PDB Queries

exports.getCathDetail = async (req, res) => {
  try {
    const id = parseInt(req.params.ID);
    let cath = await cathInfo_pdb.findOne({
      where: {
        atomic_structure_id: id
      }
    });
    if (!cath) {
      console.log("Cath information " + id + " not found.");
      res.status(400).send({
        message: "Cath information " + id + " not found."
      });
    } else {
      res.status(200).json(cath);
    }
  } catch (error) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

exports.getResultsPDB = async (req, res) => {
  try {
    const id = parseInt(req.params.ID);
    const rep = parseInt(req.params.rep);
    const type = parseInt(req.params.type);
    const cath = parseInt(req.params.cath);
    const len = parseInt(req.params.len);
    const top = parseInt(req.params.top);
    let query_results = await getBiomoleculesPDB(id, rep, type, cath, len, top);
    console.log(query_results);
    res.status(200).json(query_results);
  } catch (error) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

//#endregion

//#region Auxiliar functions

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

async function getBiomolecules(
  emd_id_p,
  type_descriptor,
  top_can,
  is_volume_filter_on,
  min_res,
  max_res
) {
  try {
    // Query to DB
    let biomolecules = await sequelize.query(
      "SELECT * FROM top_distance(:emd_id, :type_des, :top)",
      {
        replacements: {
          emd_id: emd_id_p,
          type_des: type_descriptor,
          top: top_can
        }
      }
    );
    // final result array
    let resultArray = [];
    // Clasification Process
    for (const biomoleculeItem of biomolecules[0]) {
      let bioInfo = await searchByID(biomoleculeItem["emd_id"]);
      let distance = biomoleculeItem["distance"].toString();
      if (is_volume_filter_on == "true") {
        let min_active = min_res != "0" ? true : false;
        let max_active = max_res != "0" ? true : false;
        if (min_active && !max_active && bioInfo.resolution >= min_res) {
          resultArray.push({
            biomolecule: bioInfo,
            euc_distance: distance.substring(0, 5)
          });
        } else if (!min_active && max_active && bioInfo.resolution <= max_res) {
          resultArray.push({
            biomolecule: bioInfo,
            euc_distance: distance.substring(0, 5)
          });
        } else if (
          min_active &&
          max_active &&
          bioInfo.resolution >= min_res &&
          bioInfo.resolution <= max_res
        ) {
          resultArray.push({
            biomolecule: bioInfo,
            euc_distance: distance.substring(0, 5)
          });
        }
      } else {
        resultArray.push({
          biomolecule: bioInfo,
          euc_distance: distance.substring(0, 5)
        });
      }
    }
    return resultArray;
  } catch (err) {
    return err;
  }
}

async function getBiomoleculesPDB(id_p, rep_p, db_p, cath_p, len_p, top_p) {
  try {
    // Query to DB
    console.log(id_p, rep_p, db_p, cath_p, len_p, top_p);
    let biomolecules = await sequelize.query(
      "SELECT * FROM atomic_query(:pdb_id, :rep, :db, :cath, :len, :top)",
      {
        replacements: {
          pdb_id: id_p,
          rep: rep_p,
          db: db_p,
          cath: cath_p,
          len: len_p,
          top: top_p
        }
      }
    );
    return biomolecules;
  } catch (err) {
    return err;
  }
}

function giveResponse(id, biomolecule, res) {
  if (!biomolecule) {
    console.log("Biomolecule " + id + " not found.");
    return res.status(400).send({
      message: "Biomolecule " + id + " not found."
    });
  } else {
    return res.status(200).json(biomolecule);
  }
}

//#endregion
