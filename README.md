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
git clone https://github.com/estherspiteri/Masters-project/tree/docker_compose
git checkout docker_compose 
cd vusvista
```
### 2. Build and Start the Application

Run the following command to start all services:

```bash
docker-compose up --build
```
<b>PLEASE NOTE: This step takes a while to complete since it is setting up the entire system.</b>

### 3. Access the Application

Go to [http://localhost:3001](http://localhost:3001) on your web browser.

You will be automatically redirected to the login page and in order to login you need to first register.
Since registration should be handled by authorised individuals, the registration view is inaccessible through the user interface.
<b>You need to go to [http://localhost:3001/register](http://localhost:3001/register) in order to register.</b>


- **Frontend**: [http://localhost:3001](http://localhost:3001)
- **Backend API**: [http://localhost:5001](http://localhost:5001)
- **PostgreSQL**: Runs internally on port `5432`.

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


## File Structure
Main directories and files within this repository
```plaintext
├── README.md
├── docker-compose.yml
├── .env
└── services
    ├── client
    │   ├── Dockerfile
    │   ├── README.md
    │   ├── custom.d.ts
    │   ├── nginx.conf
    │   ├── package-lock.json
    │   ├── package.json
    │   ├── public
    │   │   ├── assets
    │   │   │   └── template.xlsx
    │   │   ├── favicon.ico
    │   │   ├── index.html
    │   │   ├── manifest.json
    │   │   └── robots.txt
    │   ├── src
    │   │   ├── App.module.scss
    │   │   ├── App.tsx
    │   │   ├── app-context.tsx
    │   │   ├── assets
    │   │   │   └── logo.png
    │   │   ├── atoms
    │   │   │   ├── banner
    │   │   │   │   ├── banner.module.scss
    │   │   │   │   └── banner.tsx
    │   │   │   ├── button
    │   │   │   │   ├── button.module.scss
    │   │   │   │   └── button.tsx
    │   │   │   ├── calendar-display
    │   │   │   │   ├── calendar-display.scss
    │   │   │   │   └── calendar-display.tsx
    │   │   │   ├── dropdown
    │   │   │   │   ├── dropdown.module.scss
    │   │   │   │   └── dropdown.tsx
    │   │   │   ├── icons
    │   │   │   │   ├── icon.module.scss
    │   │   │   │   └── icon.tsx
    │   │   │   ├── loader
    │   │   │   │   ├── loader.module.scss
    │   │   │   │   └── loader.tsx
    │   │   │   ├── modal
    │   │   │   │   ├── modal.module.scss
    │   │   │   │   └── modal.tsx
    │   │   │   ├── text
    │   │   │   │   ├── text.module.scss
    │   │   │   │   └── text.tsx
    │   │   │   └── text-area
    │   │   │       ├── text-area.module.scss
    │   │   │       └── text-area.tsx
    │   │   ├── components
    │   │   │   ├── error-page
    │   │   │   │   ├── error-page.module.scss
    │   │   │   │   ├── error-page.tsx
    │   │   │   │   └── error.svg
    │   │   │   ├── header
    │   │   │   │   ├── header.module.scss
    │   │   │   │   └── header.tsx
    │   │   │   ├── home-page
    │   │   │   │   ├── home-page.module.scss
    │   │   │   │   └── home-page.tsx
    │   │   │   ├── login-page
    │   │   │   │   ├── login-page.module.scss
    │   │   │   │   └── login-page.tsx
    │   │   │   ├── profile-page
    │   │   │   │   ├── profile-page.module.scss
    │   │   │   │   └── profile-page.tsx
    │   │   │   ├── publication-view-page
    │   │   │   │   ├── publication-preview
    │   │   │   │   │   ├── publication-preview.module.scss
    │   │   │   │   │   └── publication-preview.tsx
    │   │   │   │   ├── publication-view-page.module.scss
    │   │   │   │   └── publication-view-page.tsx
    │   │   │   ├── register-page
    │   │   │   │   ├── register-page.module.scss
    │   │   │   │   └── register-page.tsx
    │   │   │   ├── review-history-page
    │   │   │   │   ├── review
    │   │   │   │   │   ├── review.module.scss
    │   │   │   │   │   └── review.tsx
    │   │   │   │   ├── review-history-page.module.scss
    │   │   │   │   └── review-history-page.tsx
    │   │   │   ├── review-page
    │   │   │   │   ├── review-page.module.scss
    │   │   │   │   └── review-page.tsx
    │   │   │   ├── sample-page
    │   │   │   │   ├── acmg-rule-info
    │   │   │   │   │   ├── acmg-rule-info.module.scss
    │   │   │   │   │   └── acmg-rule-info.tsx
    │   │   │   │   ├── acmg-rules-edit
    │   │   │   │   │   ├── acmg-rules-edit.module.scss
    │   │   │   │   │   └── acmg-rules-edit.tsx
    │   │   │   │   ├── phenotype-selection
    │   │   │   │   │   └── phenotype-selection.tsx
    │   │   │   │   ├── sample-info
    │   │   │   │   │   ├── sample-info.module.scss
    │   │   │   │   │   └── sample-info.tsx
    │   │   │   │   ├── sample-page.module.scss
    │   │   │   │   └── sample-page.tsx
    │   │   │   ├── shared
    │   │   │   │   ├── add-publications
    │   │   │   │   │   ├── add-publications.module.scss
    │   │   │   │   │   └── add-publications.tsx
    │   │   │   │   ├── sample-phenotype-selection
    │   │   │   │   │   ├── sample-phenotype-selection.module.scss
    │   │   │   │   │   └── sample-phenotype-selection.tsx
    │   │   │   │   ├── table
    │   │   │   │   │   ├── debounced-input
    │   │   │   │   │   │   ├── debounced-input.module.scss
    │   │   │   │   │   │   └── debounced-input.tsx
    │   │   │   │   │   └── filter
    │   │   │   │   │       ├── filter.module.scss
    │   │   │   │   │       └── filter.tsx
    │   │   │   │   └── variant-summary
    │   │   │   │       ├── variant-summary.module.scss
    │   │   │   │       └── variant-summary.tsx
    │   │   │   ├── view-all-vus-page
    │   │   │   │   ├── view-all-vus-page.module.scss
    │   │   │   │   ├── view-all-vus-page.tsx
    │   │   │   │   ├── view-vus
    │   │   │   │   │   ├── view-vus.module.scss
    │   │   │   │   │   └── view-vus.tsx
    │   │   │   │   └── vus-table
    │   │   │   │       ├── vus-table.module.scss
    │   │   │   │       └── vus-table.tsx
    │   │   │   ├── view-samples-page
    │   │   │   │   ├── sample-table
    │   │   │   │   │   ├── sample-table.module.scss
    │   │   │   │   │   └── sample-table.tsx
    │   │   │   │   ├── view-sample
    │   │   │   │   │   ├── view-sample.module.scss
    │   │   │   │   │   └── view-sample.tsx
    │   │   │   │   ├── view-samples-page.module.scss
    │   │   │   │   └── view-samples-page.tsx
    │   │   │   ├── vus-file-upload-page
    │   │   │   │   ├── upload.gif
    │   │   │   │   ├── vus-file-upload-page.module.scss
    │   │   │   │   └── vus-file-upload-page.tsx
    │   │   │   ├── vus-page
    │   │   │   │   ├── vus-info
    │   │   │   │   │   ├── vus-info.module.scss
    │   │   │   │   │   └── vus-info.tsx
    │   │   │   │   ├── vus-page.module.scss
    │   │   │   │   └── vus-page.tsx
    │   │   │   └── vus-upload-page
    │   │   │       ├── vus-upload-field
    │   │   │       │   ├── vus-upload-field.module.scss
    │   │   │       │   └── vus-upload-field.tsx
    │   │   │       ├── vus-upload-page.module.scss
    │   │   │       └── vus-upload-page.tsx
    │   │   ├── helpers
    │   │   │   ├── date-helper.tsx
    │   │   │   └── open-links.tsx
    │   │   ├── index.tsx
    │   │   ├── models
    │   │   │   ├── acmg-rule.model.tsx
    │   │   │   ├── classification-review.model.tsx
    │   │   │   ├── clinvar-updates.model.tsx
    │   │   │   ├── gene.model.tsx
    │   │   │   ├── homepage.model.tsx
    │   │   │   ├── load-review.model.tsx
    │   │   │   ├── phenotype.model.tsx
    │   │   │   ├── profile.model.tsx
    │   │   │   ├── publication-view.model.tsx
    │   │   │   ├── sample-to-add-info.model.ts
    │   │   │   ├── updated-external-ref-data.model.tsx
    │   │   │   ├── variant-publication-updates.tsx
    │   │   │   ├── variant-to-add-info.model.ts
    │   │   │   ├── view-samples.model.tsx
    │   │   │   ├── view-vus.model.tsx
    │   │   │   ├── vus-file-upload.model.tsx
    │   │   │   ├── vus-summary.model.tsx
    │   │   │   └── vus-upload.model.tsx
    │   │   ├── services
    │   │   │   ├── api
    │   │   │   │   └── api.service.ts
    │   │   │   ├── auth
    │   │   │   │   ├── auth.dto.ts
    │   │   │   │   └── auth.service.ts
    │   │   │   ├── homepage
    │   │   │   │   ├── homepage.dto.ts
    │   │   │   │   └── homepage.service.ts
    │   │   │   ├── profile
    │   │   │   │   ├── profile.dto.ts
    │   │   │   │   └── profile.service.ts
    │   │   │   ├── publication
    │   │   │   │   ├── publication.dto.ts
    │   │   │   │   └── publication.service.ts
    │   │   │   ├── review
    │   │   │   │   ├── review.dto.ts
    │   │   │   │   └── review.service.ts
    │   │   │   ├── sample
    │   │   │   │   ├── sample.dto.ts
    │   │   │   │   └── sample.service.ts
    │   │   │   └── vus
    │   │   │       ├── vus.dto.ts
    │   │   │       └── vus.service.ts
    │   │   ├── setupTests.js
    │   │   └── wrappers
    │   │       ├── homepage-wrapper.tsx
    │   │       ├── profile-page-wrapper.tsx
    │   │       ├── publication-phenotype-view-page-wrapper.tsx
    │   │       ├── publication-view-page-wrapper.tsx
    │   │       ├── review-history-page-wrapper.tsx
    │   │       ├── review-page-wrapper.tsx
    │   │       ├── sample-page-wrapper.tsx
    │   │       ├── vus-page-wrapper.tsx
    │   │       └── vus-upload-page-wrapper.tsx
    │   ├── tsconfig.json
    │   ├── typings.d.ts
    │   └── vercel.json
    ├── db
    │   ├── Dockerfile
    │   └── init.sql
    └── flask_server
        ├── Dockerfile
        ├── .env
        ├── README.md
        ├── app.log
        ├── entrypoint.sh
        ├── environment.yml
        ├── logs
        │   ├── app.log
        │   ├── app.log.1
        │   ├── app.log.10
        │   ├── app.log.2
        │   ├── app.log.3
        │   ├── app.log.4
        │   ├── app.log.5
        │   ├── app.log.6
        │   ├── app.log.7
        │   ├── app.log.8
        │   └── app.log.9
        ├── main.py
        ├── requirements.txt
        ├── server
        │   ├── __init__.py
        │   ├── config.py
        │   ├── db_setup
        │   │   ├── Homo_sapiens.GRCh37.87.gtf
        │   │   ├── Homo_sapiens.GRCh37.87.gtf.gz
        │   │   └── populate_gene_annotations_table.py
        │   ├── error_handlers.py
        │   ├── helpers
        │   │   ├── data_helper.py
        │   │   └── db_access_helper.py
        │   ├── models.py
        │   ├── responses
        │   │   └── internal_response.py
        │   ├── services
        │   │   ├── acmg_service.py
        │   │   ├── auth_service.py
        │   │   ├── clinvar_service.py
        │   │   ├── consequence_service.py
        │   │   ├── dbsnp_service.py
        │   │   ├── entrez_service.py
        │   │   ├── litvar_service.py
        │   │   ├── phenotype_service.py
        │   │   ├── publications_service.py
        │   │   ├── review_service.py
        │   │   ├── samples_service.py
        │   │   ├── variants_samples_service.py
        │   │   ├── view_vus_service.py
        │   │   └── vus_preprocess_service.py
        │   ├── templates
        │   │   ├── clinvar_update_email_template.html
        │   │   └── publication_update_email_template.html
        │   └── views
        │       ├── auth_views.py
        │       ├── homepage_views.py
        │       ├── profile_views.py
        │       ├── publication_views.py
        │       ├── review_views.py
        │       ├── sample_views.py
        │       └── vus_views.py
        └── variants.vcf
```

---


---

## Acknowledgments

- [Flask documentation](https://flask.palletsprojects.com/)
- [React documentation](https://reactjs.org/)
- [PostgreSQL documentation](https://www.postgresql.org/docs/)
