const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;

const descriptor = sequelize.define(
  "descriptor",
  {
    emd_entry_id: {
      type: Sequelize.INTEGER,
      primaryKey: true
    },
    type_descriptor_id: {
      type: Sequelize.INTEGER
    },
    numbers: {
      type: Sequelize.JSON
    }
  },
  {
    timestamps: false,
    freezeTableName: true,
    tableName: "descriptor"
  }
);

module.exports = descriptor;
