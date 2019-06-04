const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;

const update = sequelize.define(
  "update",
  {
    last_update: {
      type: Sequelize.DATE
    }
  },
  {
    timestamps: false,
    freezeTableName: true,
    tableName: "update"
  }
);

module.exports = update;
