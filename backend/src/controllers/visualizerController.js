const biomolecule_emd = require("../models/biomoleculeModelEMD");
const sequelize = require("../database").sequelize;
// const NGL = require("ngl/dist/ngl.esm.js")

exports.searchByID = async (req, res) => {
    try {
        const emdId = req.params.emdId;
        res.status(200).json(test_7);
      } catch (err) {
        res.status(500).send({
          message: "Backend error"
        });
      }

}

/* 
    Functions for each process 
*/



/* 
    Test variables
*/

// test arrays and colors

colors = ["white", "blue", "green", "red", "yellow", "orange", "pink", "brown",
"gray", "purple", "cyan", "silver", "bronze", "gold"];

// schemeId_1 = NGL.ColormakerRegistry.addSelectionScheme([["red", "7 or 10"], ["white", "*"]]);
// schemeId_2 = NGL.ColormakerRegistry.addSelectionScheme([["green", "7 or 10"], ["white", "*"]]);
// schemeId_3 = NGL.ColormakerRegistry.addSelectionScheme([["blue", "12 or 15"], ["white", "*"]]);
// schemeId_4 = NGL.ColormakerRegistry.addSelectionScheme([["yellow", "16 or 19"], ["white", "*"]]);
// schemeId_5 = NGL.ColormakerRegistry.addSelectionScheme([["orange", "20 or 25"], ["white", "*"]]);
// schemeId_6 = NGL.ColormakerRegistry.addSelectionScheme([["pink", "26 or 29"], ["white", "*"]]);
// schemeId_7 = NGL.ColormakerRegistry.addSelectionScheme([["purple", "31 or 35"], ["white", "*"]]);
// schemeId_8 = NGL.ColormakerRegistry.addSelectionScheme([["cyan", "36 or 41"], ["white", "*"]]);
// schemeId_9 = NGL.ColormakerRegistry.addSelectionScheme([["brown", "42 or 46"], ["white", "*"]]);

// test_1 = {0: {"path": "../../../assets/example-files/1crn/1crn.pdb", file_type: 0, color: "white", name: "1crn" }};

// test_2 = {  0:{"path":"../../../assets/example-files/1crn/1crn-1.pdb", file_type: 0, color: this.schemeId_1, name: "1crn-1"  },
//             1:{"path":"../../../assets/example-files/1crn/1crn-2.pdb", file_type: 0, color: this.schemeId_2, name: "1crn-2"  },
//             2:{"path":"../../../assets/example-files/1crn/1crn-3.pdb", file_type: 0, color: this.schemeId_3, name: "1crn-3"  },
//             3:{"path":"../../../assets/example-files/1crn/1crn-4.pdb", file_type: 0, color: this.schemeId_4, name: "1crn-4"  },
//             4:{"path":"../../../assets/example-files/1crn/1crn-5.pdb", file_type: 0, color: this.schemeId_5, name: "1crn-5"  },
//             5:{"path":"../../../assets/example-files/1crn/1crn-6.pdb", file_type: 0, color: this.schemeId_6, name: "1crn-6"  },
//             6:{"path":"../../../assets/example-files/1crn/1crn-7.pdb", file_type: 0, color: this.schemeId_7, name: "1crn-7"  },
//             7:{"path":"../../../assets/example-files/1crn/1crn-8.pdb", file_type: 0, color: this.schemeId_8, name: "1crn-8"  },
//             8:{"path":"../../../assets/example-files/1crn/1crn-9.pdb", file_type: 0, color: this.schemeId_9, name: "1crn-9"  }};

// test_3 = {0: {"path": "../../../assets/example-files/emd_1884.map", file_type: 1, color: "white", name: "emd_1884"  }};

// test_4 = {0: {"path": "../../../assets/example-files/175d/175d.mrc", file_type: 1, color: "white", name: "175d"  }};

// test_5 = {  0:{"path":"../../../assets/example-files/175d/175d_A.mrc", file_type: 1, color: "white", name: "175d_A"  },
//             1:{"path":"../../../assets/example-files/175d/175d_B.mrc", file_type: 1, color: "white", name: "175d_B"  }};

// test_6 = {0: {"path": "../../../assets/example-files/1crn_mrc/1crn.mrc", file_type: 1, color: "white", name: "1crn"  }};

test_7 = {  0:{"path":"../../../assets/example-files/1crn_mrc/1crn_A.mrc", file_type: 1, color: colors[0], name: "1crn_A" },
            1:{"path":"../../../assets/example-files/1crn_mrc/1crn_B.mrc", file_type: 1, color: colors[1], name: "1crn_B"  },
            2:{"path":"../../../assets/example-files/1crn_mrc/1crn_C.mrc", file_type: 1, color: colors[2], name: "1crn_C"  },
            3:{"path":"../../../assets/example-files/1crn_mrc/1crn_D.mrc", file_type: 1, color: colors[3], name: "1crn_D"  },
            4:{"path":"../../../assets/example-files/1crn_mrc/1crn_E.mrc", file_type: 1, color: colors[4], name: "1crn_E"  },
            5:{"path":"../../../assets/example-files/1crn_mrc/1crn_F.mrc", file_type: 1, color: colors[5], name: "1crn_F"  },
            6:{"path":"../../../assets/example-files/1crn_mrc/1crn_G.mrc", file_type: 1, color: colors[6], name: "1crn_G"  },
            7:{"path":"../../../assets/example-files/1crn_mrc/1crn_H.mrc", file_type: 1, color: colors[7], name: "1crn_H"  },
            8:{"path":"../../../assets/example-files/1crn_mrc/1crn_I.mrc", file_type: 1, color: colors[8], name: "1crn_I"  },
            9:{"path":"../../../assets/example-files/1crn_mrc/1crn_J.mrc", file_type: 1, color: colors[9], name: "1crn_J"  },
            10:{"path":"../../../assets/example-files/1crn_mrc/1crn_K.mrc", file_type: 1, color: colors[10], name: "1crn_K"  },
            11:{"path":"../../../assets/example-files/1crn_mrc/1crn_L.mrc", file_type: 1, color: colors[11], name: "1crn_L"  },
            12:{"path":"../../../assets/example-files/1crn_mrc/1crn_M.mrc", file_type: 1, color: colors[12], name: "1crn_M"  },
            13:{"path":"../../../assets/example-files/1crn_mrc/1crn_N.mrc", file_type: 1, color: colors[13], name: "1crn_N"  }};
