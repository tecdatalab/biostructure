const args = process.argv.slice(2);
const Sequelize = require("sequelize");
const sequelize = new Sequelize(
  args[4], //db name
  args[0], //username
  args[1], //password
  {
    dialect: "postgres",
    host: args[2],
    port: args[3],
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
