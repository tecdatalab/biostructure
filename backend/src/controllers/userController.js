const googleAuth = require("../security/googleAuth");
const jwt = require("../security/jwt");
const user = require("../models/userModel");
const userRole = require("../models/userRoleModel");

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
    res.status(500).send({ message: "Internal server error" });
    return console.error(error);
  }
};

exports.verifyUserToken = async (req, res, next) => {
  console.log("Verify");
  const token = req.header("authorization");
  jwt.verify(token, (err, decodedToken) => {
    if (err) {
      res.status(401).send({ message: "user is not authenticated" });
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
            res.status(401).send({ message: "user is not authenticated" });
          }
        })
        .catch(error => {
          res.status(403).send({
            message: "user is not authorized or does not exists"
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
      .send({ message: "user is not authorized or does not exists" });
  }
};

exports.changeUserRole = (req, res) => {
  const userId = req.body.userId;
  const newUserRole = req.body.role;
  user
    .update({ role: newUserRole }, { where: { id: userId }, returning: true })
    .then(result => {
      if (result[0] != 0) {
        res.status(200).send(result[1][0]);
      } else {
        throw new Error("Could not assign admin role to the requested user");
      }
    })
    .catch(error => {
      res.status(400).send({ message: error });
      return console.error(error);
    });
};

exports.getUsersRoles = async (req, res) => {
  try {
    let userRoles = await userRole.findAll({});
    res.status(200).json(userRoles);
  } catch (err) {
    res.status(500).send({
      message: "Internal server error"
    });
  }
};

exports.getUsers = async (req, res) => {
  try {
    let users = await user.findAll({});
    res.status(200).json(users);
  } catch (err) {
    res.status(500).send({
      message: "Internal server error"
    });
  }
};

exports.isUserAdmin = (req, res) => {
  const reqUser = req.authenticatedUser;
  if (reqUser.role == 2) {
    res.status(200).send(true);
  } else {
    res.status(200).send(false);
  }
};

exports.isUserLoggedIn = (req, res, next) => {
  if (req.header("authorization")) {
    exports.verifyUserToken(req, res, next);
  } else {
    req.authenticatedUser = { id: null };
    next();
  }
};
