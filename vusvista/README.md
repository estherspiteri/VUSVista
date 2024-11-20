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
- [Docker Compose](https://docs.docker.com/compose/)

---

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/estherspiteri/Masters-project/tree/docker_compose
cd vusvista
```

### 2. Create Environment File

In the root directory, create a `.env` file based on `.env.example` but fill in the correct values.

### 3. Build and Start the Application

Run the following command to start all services:

```bash
docker-compose up --build
```

### 4. Access the Application

Go to [http://localhost:3000](http://localhost:3000) on your web browser.

- **Frontend**: [http://localhost:3000](http://localhost:3000)
- **Backend API**: [http://localhost:5000](http://localhost:5000)
- **PostgreSQL**: Runs internally on port `5432`.

---

## File Structure

```plaintext
├── backend
│   ├── app.py               # Flask app entry point
│   ├── requirements.txt     # Backend dependencies
│   └── Dockerfile           # Backend Docker configuration
├── frontend
│   ├── src/
│   ├── package.json         # Frontend dependencies
│   └── Dockerfile           # Frontend Docker configuration
├── docker-compose.yml       # Multi-service orchestration
└── README.md                # Project documentation
```

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

You can connect to the PostgreSQL database using tools like `psql` or a GUI tool such as **DBeaver** or **pgAdmin**. Use the credentials provided in the `docker-compose.yml` file.

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

- **Port Conflicts**: Ensure no other services are running on ports `3000` (React), `5000` (Flask), or `5432` (PostgreSQL).
- **Database Issues**: Check the `DATABASE_URL` in the `.env` file.

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

---

## Acknowledgments

- [Flask documentation](https://flask.palletsprojects.com/)
- [React documentation](https://reactjs.org/)
- [PostgreSQL documentation](https://www.postgresql.org/docs/)
