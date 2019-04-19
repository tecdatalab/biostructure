const jwt = require("jsonwebtoken");

const secret = require("./credentials.json").jwt.secret; // TO DO: create an ENV variable to contain the secret

exports.verify = token => {
  try {
    return jwt.verify(token, secret);
  } catch (e) {
    throw new Error("jwt token not verified");
  }
};

exports.generateToken = user => {
  expDate = new Date();
  expDate.setDate(expDate.getDate() + 8); //TO DO: define the expiration time in a config file
  return jwt.sign(
    {
      sub: user.id,
      exp: expDate.getTime() // expiration date given in milliseconds
    },
    secret
  );
};
