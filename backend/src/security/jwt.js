const jwt = require("jsonwebtoken");

const secret = require("./credentials.json").jwt.secret; // TO DO: create an ENV variable to contain the secret

exports.verify = (token, errFunc) => {
  return jwt.verify(token, secret, errFunc);
};

exports.generateToken = user => {
  expDate = new Date();
  expDate.setDate(expDate.getDate() + 8); //TO DO: define the expiration time in a config file
  return jwt.sign(
    {
      user: user,
      exp: expDate.getTime() // expiration date given in milliseconds
    },
    secret
  );
};
