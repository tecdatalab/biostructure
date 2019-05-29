/*
 * Author: José Antonio Salas-Bonilla, Julio Viquez-Murillo.
 * Purpose: Creation of tables for the management of users.
 *
 */

CREATE TABLE user_role(
    id SERIAL PRIMARY KEY,
    role TEXT NOT NULL
);

COMMENT ON TABLE user_role IS 'Stores the different types of user roles ';

CREATE TABLE "user"(
    id TEXT PRIMARY KEY,
    name TEXT NOT NULL,
    email TEXT NOT NULL UNIQUE,
    role INT REFERENCES user_role(id) DEFAULT 1 NOT NULL
);

COMMENT ON TABLE "user" IS 'Stores the user information';

