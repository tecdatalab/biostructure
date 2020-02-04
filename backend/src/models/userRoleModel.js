const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;

const userRole = sequelize.define(
  "user_role",
  {
    id: {
      type: Sequelize.TEXT,
      primaryKey: true
    },
    role: {
      type: Sequelize.INTEGER
    }
  },
  {
    timestamps: false,
    freezeTableName: true,
    tableName: "user_role"
  }
);

module.exports = userRole;
