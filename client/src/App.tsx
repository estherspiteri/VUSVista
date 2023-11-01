import React from "react";
import styles from "./App.module.scss";
import PublicationSearch from "./components/publication-search/publication-search";
import { publicationService } from "./services/publication/publication.service";

type AppProps = {};

const App: React.FunctionComponent<AppProps> = () => {
  return (
    <div className={styles["container"]}>
      <PublicationSearch publicationService={publicationService} />
    </div>
  );
};

export default App;
