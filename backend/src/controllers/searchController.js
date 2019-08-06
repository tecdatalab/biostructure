const biomolecule = require("../models/biomoleculeModel");
const search_history = require("../models/searchHistoryModel");
const Op = require("../database").Op;
const sequelize = require("../database").sequelize;

exports.searchByID = async (req, res) => {
  try {
    const emdbid = req.params.emdbID;
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

    // let query_results = await getBiomolecules(minRes, maxRes);
    let query_results = await getBiomolecules(12, 1, 10);

    console.log("\n");
    console.log(query_results);
    console.log("\n")

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

/**
 * Function 'getBiomolecules' Updated to work with the SP of Luis in the DB
 * Updated by: Esteban Sanabria
 * Date: 24/07/2019
 */

// OLD VERSION
// async function getBiomolecules(minRes, maxRes) {
//   try {
//     let biomolecules = await biomolecule.findAll({
//       where: {
//         volume: {
//           [Op.gte]: minRes,
//           [Op.lte]: maxRes
//         }
//       }
//     });
//     let resultArray = [];
//     biomolecules.forEach(biomoleculeItem => {
//       resultArray.push({
//         biomolecule: biomoleculeItem.get({
//           plain: true
//         }),
//         euc_distance: 5
//       });
//     });
//     return resultArray;
//   } catch (err) {
//     return err;
//   }
// }

async function getBiomolecules(emd_id_p, type_descriptor, top_can) {
  try {
    // Query to DB
    let biomolecules = await sequelize.query(
                  'SELECT * FROM top_distance(:emd_id, :type_des, :top)',
                  {replacements: { emd_id: emd_id_p, 
                                   type_des: type_descriptor, 
                                   top: top_can }
                                  });
    // final result array
    let resultArray = [];
    // Clasification Process
    biomolecules[0].forEach(biomoleculeItem => {
      resultArray.push({
                biomolecule: searchByID(biomoleculeItem["emd_id"]),
                euc_distance: biomoleculeItem["distance"]
              });
    });     
    return resultArray;
  } catch (err) {
    return err;
  }
}

async function searchByID(emdbid){
  try {
    await biomolecule.findOne({
      where: {
        id: parseInt(emdbid)
      }
    }).then(biomolecules => {
        if (!biomolecules) {
          console.log("Biomolecule " + emdbid + " not found.");
          return -1;
        } else {
          return biomolecules.get({plain: true});
        }
      }
    );
  } catch (err) {
    return -1;
  }
}