const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;

const biomolecule = sequelize.define(
    "pdb_entry",
    {
        id: {
            type: Sequelize.INTEGER,
            primaryKey : true
        },
        id_code: {
            type: Sequelize.STRING
        },
        numbers_descriptor: {
            type: Sequelize.JSON
        },
        sequence: {
            type: Sequelize.STRING
        },
        png_img_3d: {
            type: Sequelize.STRING
        },
        gif_img_3d: {
            type: Sequelize.STRING
        }
    },
    {
        timestamps: false,
        freezeTableName: true,
        tableName: "atomic_structure"
    }
);

module.exports = biomolecule;
