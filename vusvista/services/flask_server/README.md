# VUSVista - Flask Project Setup

This project is a web application built with Flask (backend), PostgreSQL (database), and React (frontend). Below are instructions to run the Flask backend locally for development purposes.

---

## Prerequisites

Ensure the following tools are installed:

- [Python 3](https://www.python.org/downloads/)
- [pip](https://pip.pypa.io/en/stable/)
- [Virtualenv](https://virtualenv.pypa.io/en/stable/)

---

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/estherspiteri/Masters-project/tree/docker_compose
cd vusvista
```

### 2. Navigate to the Backend Directory

```bash
cd backend
```

### 3. Create a Virtual Environment

Run the following command to create a virtual environment:

```bash
python -m venv venv
```

### 4. Activate the Virtual Environment

- **Linux/macOS**:
  ```bash
  source venv/bin/activate
  ```
- **Windows**:
  ```bash
  venv\Scripts\activate
  ```

### 5. Install Dependencies

Install the required Python packages:

```bash
pip install -r requirements.txt
```

### 6. Create Environment File

In the `backend` directory, create a `.env` file based on `.env.example` but fill in the correct values.

### 7. Ensure the Database is Running
Checkout `README.md` file in `vusvista/services/db`

### 8. Start the Flask Application

Run the following command to start the Flask application locally:

```bash
flask run
```

The application will start in development mode, and you can access the API at [http://localhost:5000](http://localhost:5000).

---

## Development Workflow

- **Running the Backend Locally**:
  ```bash
  cd backend
  flask run
  ```

- **Install New Dependencies**:
  If you need additional Python packages, install them and save to `requirements.txt`:
  ```bash
  pip install <package-name>
  pip freeze > requirements.txt
  ```

---

## Useful Commands

### Stop the Flask Development Server:

Use `CTRL + C` in the terminal to stop the development server.

---

## Troubleshooting

- **Port Conflicts**: Ensure no other services are running on port `5000` (Flask).
- **Missing Dependencies**: Run `pip install -r requirements.txt` if you encounter errors related to missing packages.

---

## Acknowledgments

- [Flask documentation](https://flask.palletsprojects.com/)

