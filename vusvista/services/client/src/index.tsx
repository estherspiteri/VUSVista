import React from "react";
import ReactDOM from "react-dom";
import App from "./App.tsx";
import { BrowserRouter } from "react-router-dom";
import { AppProvider } from "./app-context.tsx";

ReactDOM.render(
  <React.StrictMode>
    <BrowserRouter>
      <AppProvider>
        <App />
      </AppProvider>
    </BrowserRouter>
  </React.StrictMode>,
  document.getElementById("root")
);
