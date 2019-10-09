const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;

const cath = sequelize.define(
    "cathDetail",
    {
      cath_domain_name: {
          type: Sequelize.TEXT
      },
      class_number: {
          type: Sequelize.INTEGER
      },
      arquitecture_number: {
          type: Sequelize.INTEGER
      },
      topology_number: {
          type: Sequelize.INTEGER
      },
      homologous_superfamily_number: {
          type: Sequelize.INTEGER
      },
      atomic_structure_id: {
          type: Sequelize.INTEGER
      }
    },
    {
      timestamps: false,
      freezeTableName: true,
      tableName: "cath_x_atomic_structure"
    }
);

module.exports = cath;