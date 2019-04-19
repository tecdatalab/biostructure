const googleAuth = require("../security/googleAuth");
const jwt = require("../security/jwt");

const getCredentials = login => {
  return googleAuth.getGoogleUser(login.code).then(googleUser => {
    const content = {
      token: jwt.generateToken(googleUser),
      user: googleUser
    };
    return content;
  });
};

exports.sendAuthToken = async (req, res, next) => {
  try {
    const login = req.body;
    getCredentials
      .then(credentials => {
        res.json(credentials).end();
      })
      .catch(e => {
        console.log(e);
        throw new Error(e);
      });
  } catch (error) {
    res.sendStatus(500).end(JSON.stringify({ error: "Internal server error" }));
    return console.error(error);
  }
};
