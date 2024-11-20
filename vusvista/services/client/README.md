# VUSVista - Run React Locally

This project is a web application built with Flask (backend), PostgreSQL (database), and React (frontend). Below are instructions to run the React frontend locally for development purposes.

---

## Prerequisites

Ensure the following tools are installed:

- [Node.js and npm](https://nodejs.org/)

---

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/estherspiteri/Masters-project/tree/docker_compose
git checkout docker_compose 
cd vusvista
```

### 2. Navigate to the Frontend Directory

```bash
cd services
cd client
```

### 3. Install Dependencies

Run the following command to install all the necessary dependencies:

```bash
npm install
```
### 4. Create Environment File

In the root directory, create a `.env` file based on `.env.example` but fill in the correct values.


### 5. Start the Frontend Application

Run the following command to start the React application locally:

```bash
npm start
```

This will start the React development server, and you can view the application at [http://localhost:3000](http://localhost:3000).

---

## Development Workflow

- **Running the Frontend Locally**:
  ```bash
  cd frontend
  npm start
  ```

- **Install New Dependencies**:
  If you need additional npm packages, install them and save to `package.json`:
  ```bash
  npm install <package-name>
  ```

---

# Getting Started with Create React App

This project was bootstrapped with [Create React App](https://github.com/facebook/create-react-app).

## Available Scripts

In the project directory, you can run:

### `npm start`

Runs the app in the development mode.\
Open [http://localhost:3000](http://localhost:3000) to view it in your browser.

The page will reload when you make changes.\
You may also see any lint errors in the console.

### `npm test`

Launches the test runner in the interactive watch mode.\
See the section about [running tests](https://facebook.github.io/create-react-app/docs/running-tests) for more information.

### `npm run build`

Builds the app for production to the `build` folder.\
It correctly bundles React in production mode and optimizes the build for the best performance.

The build is minified and the filenames include the hashes.\
Your app is ready to be deployed!

See the section about [deployment](https://facebook.github.io/create-react-app/docs/deployment) for more information.

### `npm run eject`

**Note: this is a one-way operation. Once you `eject`, you can't go back!**

If you aren't satisfied with the build tool and configuration choices, you can `eject` at any time. This command will remove the single build dependency from your project.

Instead, it will copy all the configuration files and the transitive dependencies (webpack, Babel, ESLint, etc) right into your project so you have full control over them. All of the commands except `eject` will still work, but they will point to the copied scripts so you can tweak them. At this point you're on your own.

You don't have to ever use `eject`. The curated feature set is suitable for small and middle deployments, and you shouldn't feel obligated to use this feature. However we understand that this tool wouldn't be useful if you couldn't customize it when you are ready for it.

## Learn More

You can learn more in the [Create React App documentation](https://facebook.github.io/create-react-app/docs/getting-started).

To learn React, check out the [React documentation](https://reactjs.org/).

### Code Splitting

This section has moved here: [https://facebook.github.io/create-react-app/docs/code-splitting](https://facebook.github.io/create-react-app/docs/code-splitting)

### Analyzing the Bundle Size

This section has moved here: [https://facebook.github.io/create-react-app/docs/analyzing-the-bundle-size](https://facebook.github.io/create-react-app/docs/analyzing-the-bundle-size)

### Making a Progressive Web App

This section has moved here: [https://facebook.github.io/create-react-app/docs/making-a-progressive-web-app](https://facebook.github.io/create-react-app/docs/making-a-progressive-web-app)

### Advanced Configuration

This section has moved here: [https://facebook.github.io/create-react-app/docs/advanced-configuration](https://facebook.github.io/create-react-app/docs/advanced-configuration)

### Deployment

This section has moved here: [https://facebook.github.io/create-react-app/docs/deployment](https://facebook.github.io/create-react-app/docs/deployment)

### `npm run build` fails to minify

This section has moved here: [https://facebook.github.io/create-react-app/docs/troubleshooting#npm-run-build-fails-to-minify](https://facebook.github.io/create-react-app/docs/troubleshooting#npm-run-build-fails-to-minify)

---

## Useful Commands

### Stop the React Development Server:

Use `CTRL + C` in the terminal to stop the development server.

### Build the Frontend for Production:

```bash
npm run build
```

This command will create an optimized production build in the `build` directory.

---

## Troubleshooting

- **Port Conflicts**: Ensure no other services are running on port `3000` (React).
- **Missing Dependencies**: Run `npm install` if you encounter errors related to missing packages.


---

## Acknowledgments

- [React documentation](https://reactjs.org/)