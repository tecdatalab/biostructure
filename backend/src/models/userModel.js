const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;

const user = sequelize.define(
  "user",
  {
    id: {
      type: Sequelize.TEXT,
      primaryKey: true
    },
    name: {
      type: Sequelize.TEXT
    },
    email: {
      type: Sequelize.TEXT
    },
    role: {
      type: Sequelize.INTEGER
    }
  },
  {
    timestamps: false,
    freezeTableName: true,
    tableName: "user"
  }
);

module.exports = user;
