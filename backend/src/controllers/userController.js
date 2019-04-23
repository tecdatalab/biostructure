const googleAuth = require("../security/googleAuth");
const jwt = require("../security/jwt");
const user = require("../models/userModel");

const getGoogleUser = login => {
  return googleAuth.getGoogleUser(login.tokenId).then(googleUser => {
    return googleUser;
  });
};

exports.sendAuthToken = async (req, res) => {
  try {
    const login = req.body;
    getGoogleUser(login)
      .then(googleUser => {
        user
          .findOne({
            attributes: ["id", "name", "email"],
            where: { id: googleUser.id }
          })
          .then(dbUser => {
            if (!dbUser) {
              user
                .build({
                  id: googleUser.id,
                  name: googleUser.name,
                  email: googleUser.email
                })
                .save()
                .then(newUser => {
                  const response = {
                    token: jwt.generateToken(newUser),
                    user: newUser
                  };
                  res.json(response).end();
                })
                .catch(e => {
                  throw new Error(e);
                });
            } else {
              const response = {
                token: jwt.generateToken(dbUser),
                user: dbUser
              };
              res.json(response).end();
            }
          });
      })
      .catch(e => {
        console.log(e);
        throw new Error(e);
      });
  } catch (error) {
    res.status(500).send({ Error: "Internal server error" });
    return console.error(error);
  }
};

exports.verifyUserToken = (req, res, next) => {
  const token = req.body.token;
  jwt.verify(token, (err, decodedToken) => {
    if (err) {
      res.status(401).send({ Error: "user is not authenticated" });
    } else {
      user
        .findOne({
          where: {
            id: decodedToken.user.id
          }
        })
        .then(result => {
          if (result) {
            req.authenticatedUser = result;
            next();
          } else {
            res.status(401).send({ Error: "user is not authenticated" });
          }
        })
        .catch(error => {
          res.status(403).send({
            Error: "user is not authorized or does not exists"
          });
        });
    }
  });
};

exports.verifyAdminUser = (req, res, next) => {
  const reqUser = req.authenticatedUser;
  if (reqUser.role == 2) {
    next();
  } else {
    res
      .status(403)
      .send({ Error: "user is not authorized or does not exists" });
  }
};

exports.grantAdminRole = (req, res, next) => {
  const userId = req.body.userId;
  user
    .update({ role: 2 }, { where: { id: userId }, returning: true })
    .then(result => {
      if (result[0] != 0) {
        res.status(200).send(result[1][0]);
      } else {
        throw new Error("Could not assign admin role to the requested user");
      }
    })
    .catch(error => {
      res.status(400).send({ Error: error });
      return console.error(error);
    });
};
