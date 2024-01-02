import React, { useState } from "react";
import styles from "./App.module.scss";
import { IVus } from "./models/view-vus.model.tsx/view-vus.model";
import { Route, Routes } from "react-router-dom";
import PublicationSearch from "./components/publication-search-page/publication-search-page";
import { publicationService } from "./services/publication/publication.service";
import { vusService } from "./services/vus/vus.service";
import Header from "./components/header/header";
import ViewVusPage from "./components/view-vus-page/view-vus-page";
import VusFileUploadPage from "./components/vus-file-upload-page/vus-file-upload-page";
import ViewAllSamples from "./components/view-samples-page/view-samples-page";

type AppProps = {};
//TODO: add session cookie Id
const App: React.FunctionComponent<AppProps> = () => {
  return (
    //TODO: lazy loading
    //routing: https://hygraph.com/blog/routing-in-react
    <div className={styles.container}>
      <Header />
      <Routes>
        <Route
          path="/upload-vus"
          element={<VusFileUploadPage vusService={vusService} />}
        />
        <Route
          path="/view-samples"
          element={<ViewAllSamples vusService={vusService} />}
        />
        <Route
          path="/view-vus"
          element={<ViewVusPage vusService={vusService} />}
        />
        <Route
          path="/publication-search/*"
          element={
            <PublicationSearch publicationService={publicationService} />
          }
        />
        {/*TODO: handle no route match*/}
        {/* <Route path="*" element={<NoMatch />} /> */}
      </Routes>
    </div>
  );
};

export default App;
