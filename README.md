# VUSVista

### A Flask-PostgreSQL-React System Run with Docker Compose

This project is a web application built with Flask (backend), PostgreSQL (database), and React (frontend). The system is containerized and orchestrated using Docker Compose for ease of deployment and development.

---

## Features

- **Backend**: Flask RESTful API
- **Frontend**: React Single Page Application (SPA)
- **Database**: PostgreSQL
- **Containerization**: Docker Compose for simplified deployment

---

## Prerequisites

Ensure the following tools are installed:

- [Docker](https://www.docker.com/)
- [Docker Compose v2.x](https://docs.docker.com/compose/)
(specifically v2.x because the `depends_on` is unavailable in v3.x)
---

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/estherspiteri/VUSVista.git
cd VUSVista
```
### 2. Build and Start the Application

Run the following command to start all services:

```bash
docker-compose up --build
```
<b>PLEASE NOTE: This step takes a while to complete since it is setting up the entire system.</b>

### 3. Access the Application

- **Frontend**: [http://localhost:3001](http://localhost:3001)
- **Backend API**: [http://localhost:5001](http://localhost:5001)
- **PostgreSQL**: Runs internally on port `5432`.

Go to [http://localhost:3001](http://localhost:3001) on your web browser.

You will be automatically redirected to the login page and in order to login you need to first register.
Since registration should be handled by authorised individuals, the registration view is inaccessible through the user interface.
<b>You need to go to [http://localhost:3001/register](http://localhost:3001/register) in order to register.</b>

#### However, for the purpose of testing the system you can also login using a pre-existing user:</br>
- Email: vus.curation.system@gmail.com</br>
- Password: demoUser123

---

## Development Workflow

### Running Individual Components

#### Backend:

```bash
cd backend
flask run
```

#### Frontend:

```bash
cd frontend
npm start
```

### Accessing the Database

You can connect to the PostgreSQL database using tools like `psql` or a GUI tool such as **pgAdmin**. Use the credentials provided in the `docker-compose.yml` file. The host name/address field needs to be set to `localhost`.

---

## Useful Commands

### Stop Services:

```bash
docker-compose down
```

### View Logs:

```bash
docker-compose logs -f
```

### Rebuild Services:

```bash
docker-compose up --build
```

---

## Troubleshooting

- **Port Conflicts**: Ensure no other services are running on ports `3001` (React), `5001` (Flask), or `5432` (PostgreSQL).
- **Database Issues**: Check the `DATABASE_URL` in the `.env` file.


---

## Acknowledgments

- [Flask documentation](https://flask.palletsprojects.com/)
- [React documentation](https://reactjs.org/)
- [PostgreSQL documentation](https://www.postgresql.org/docs/)
