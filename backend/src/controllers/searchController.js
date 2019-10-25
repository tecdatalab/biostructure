const biomolecule_emd = require("../models/biomoleculeModelEMD");
const biomolecule_pdb = require("../models/biomoleculeModelPDB");
const cathInfo_pdb = require("../models/cathModel");
const search_history = require("../models/searchHistoryModel");
const Op = require("../database").Op;

//#region General Functions

exports.searchByID = async (req, res) => {
  try {
    const table = req.params.app;
    let id = req.params.ID;
    let biomolecule = null;
    if (table=="emdb"){
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
  }catch (error) {
    res.status(500).send({
      message: "Backend error"
    });
  }
}

exports.getResultsPDB = async (req, res) => {
  try {
    const c = parseInt(req.params.C);
    const a = parseInt(req.params.A);
    const t = parseInt(req.params.T);
    const h = parseInt(req.params.H);
    const id = parseInt(req.params.ID);
    let query_results = getBiomoleculesPDB(c,a,t,h,id);
    res.status(200).json(query_results);
  } catch (error) {
    res.status(500).send({
      message: "Backend error"
    });
  }
  
}

//#endregion

//#region Auxiliar functions

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
    let biomolecules = await biomolecule_emd.findAll({
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

async function getBiomoleculesPDB(c,a,t,h,id_pdb){
  try {
    let caths = await cathInfo_pdb.findAll({
      where : {
        class_numbervolume : c,
        architecture_number : a,
        topology_number : t,
        homologous_superfamily_number : h,
        atomic_structure_id: { $not : id_pdb }
      }
    });
    let biomolecules = [];
    caths.forEach(cathItem => {
      biomolecules.push({
        biomolecule: biomoleculeItem.findOne({
          where :{
            id : cathItem.atomic_structure_id
          }
      })
    });
  });
  let resultArray = [];
  for (let index = 0; index < biomolecules.length; index++) {
    resultArray.push(biomolecules[index],caths[index]);  
  }
    return resultArray;
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