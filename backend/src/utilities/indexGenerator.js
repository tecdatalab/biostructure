exports.createIndex = function() {
  return (
    new Date().valueOf().toString() + Math.floor(Math.random() * 100).toString()
  );
};
