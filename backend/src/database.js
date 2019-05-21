const Sequelize = require("sequelize");
const sequelize = new Sequelize(
  "biomolecules_db", //db name
  "postgres", //username
  "root", //password
  {
    dialect: "postgres",
    host: "localhost",
    omitNull: true,
    operatorAliases: false,
    pool: {
      max: 5,
      min: 0,
      require: 30000,
      idle: 10000
    }
  }
);

sequelize
  .authenticate()
  .then(() => console.log("Database conected"))
  .catch(err => console.log("Error " + err));

const Op = Sequelize.Op;
module.exports = {
  sequelize,
  Op
};
