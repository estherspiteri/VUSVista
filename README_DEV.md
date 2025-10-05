# VUSVista (Running components seperately i.e. not using Docker Compose)

### A Flask-PostgreSQL-React System

This project is a web application built with Flask (backend), PostgreSQL (database), and React (frontend).

---

## Features

- **Backend**: Flask RESTful API
- **Frontend**: React Single Page Application (SPA)
- **Database**: PostgreSQL

---

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/estherspiteri/VUSVista.git
cd VUSVista/vusvista
```

### 2. Set up environment variables

#### - Database
1. Create a PostgreSQL database. You can connect to the PostgreSQL database using tools like `psql` or a GUI tool such as **pgAdmin**.
2. Run the script in `vusvista/services/db/init.sql` in your database

#### - Flask Server
1. Create a copy of the file `vusvista/services/flask_server/env.example` in the `vusvista/services/flask_server` directory and name it `.env`. 
2. For the variable `SQLALCHEMY_DATABASE_URI` change the host name/address field (i.e. `db`) to `localhost`
3. Populate missing variables and update as desired

#### - React Front
1. Create a copy of the file `vusvista/services/client/env.example` in the `vusvista/services/client` directory and name it `.env`. Make sure to set `NODE_ENV` to `development`. `REACT_APP_API_URL` is not used in this scenario since the `proxy` field in the `package.json` file is used instead to connect to the backend.
2. Update variables as desired

## IMP: For both the Flask Server and the React Front, if you already have a `.env` file for the docker compose setup (`READ_ME.md`), you can create a `.env.local` file in the respective directory which overrides the variables in your `.env` file.

### 3. Build and Start the Application

Run the following command to start the Flask Server:

```bash
cd services/flask_server && flask run
```
Run the following command to start the React Front:

```bash
cd services/client && npm run start
```

### 4. Access the Application

- **Frontend**: [http://localhost:3000](http://localhost:3000)
- **Backend API**: [http://localhost:5000](http://localhost:5000)
- **PostgreSQL**: Runs internally on port `5432`.

Go to [http://localhost:3000](http://localhost:3000) on your web browser.

You will be automatically redirected to the login page and in order to login you need to first register.
Since registration should be handled by authorised individuals, the registration view is inaccessible through the user interface.
<b>You need to go to [http://localhost:3000/register](http://localhost:3000/register) in order to register.</b>

#### However, for the purpose of testing the system you can also login using a pre-existing user:</br>
- Email: vus.curation.system@gmail.com</br>
- Password: demoUser123
#### If you don't want this user to be created, please comment out or remove lines 81-83 in `vusvista/services/flask_server/server/__init__.py` before loading the system for the first time!

---

## Troubleshooting

- **Port Conflicts**: Ensure no other services are running on ports `3000` (React), `5000` (Flask), or `5432` (PostgreSQL).

---

## Acknowledgments

- [Flask documentation](https://flask.palletsprojects.com/)
- [React documentation](https://reactjs.org/)
- [PostgreSQL documentation](https://www.postgresql.org/docs/)
