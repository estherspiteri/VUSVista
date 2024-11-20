# VUSVista - PostgreSQL Database Setup

This README provides instructions for setting up the PostgreSQL database, which includes initializing the schema using the provided `init.sql` script. The database is initially set up using this script and continues to be populated when the Flask backend is started.

---

## Prerequisites

Ensure the following tools are installed:

- [PostgreSQL](https://www.postgresql.org/download/) (to run PostgreSQL locally)

---

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/estherspiteri/Masters-project/tree/docker_compose
git checkout docker_compose 
cd vusvista
```

### 2. Navigate to the PostgreSQL Directory

```bash
cd services
cd db
```

### 3. Start the PostgreSQL Service Locally

Make sure that PostgreSQL is installed and running on your machine.

### 4. Create the Database

Log in to PostgreSQL and create the database:

```bash
psql -U postgres
```

Once logged in, run the following SQL command to create the database:

```sql
CREATE DATABASE vcs_db;
```

### 5. Run the Initialization Script

Execute the `init.sql` script to set up the database schema:

```bash
psql -U postgres -d vcs_db -f init.sql
```

This will create the necessary tables and initialize the schema.

### 6. Verify Database Initialization

You can connect to the PostgreSQL database to verify that the schema was created correctly:

```bash
psql -U postgres -d vcs_db
```

Run SQL commands to verify that the tables were created successfully.

---

## Notes

- **Database Population**: The database is initially created using `init.sql`. Additional data is populated when the Flask backend is started.
- **Environment Variables**: Ensure the `.env` file contains the correct database credentials and settings for the backend to connect successfully.

---

## Useful Commands

### Start PostgreSQL Service:

Ensure that the PostgreSQL service is running. On Linux systems, you can use:

```bash
sudo service postgresql start
```

### Connect to the PostgreSQL Database:

```bash
psql -U postgres -d vcs_db
```

---

## Troubleshooting

- **Database Connection Issues**: Ensure that the PostgreSQL service is running, and check the connection settings in your `.env` file.
- **Port Conflicts**: Make sure that port `5432` is available and not used by other services.

---

## Acknowledgments

- [PostgreSQL documentation](https://www.postgresql.org/docs/)

