import React, { useState } from "react";
import styles from "./App.module.scss";
import PublicationSearch from "./components/publication-search/publication-search";
import { publicationService } from "./services/publication/publication.service";
import VusFileUpload from "./components/vus-file-upload/vus-file-upload";
import { vusService } from "./services/vus/vus.service";
import { IVus } from "./models/view-vus.model.tsx/view-vus.model";
import ViewAllVus from "./components/view-all-vus/view-all-vus";

type AppProps = {};

const App: React.FunctionComponent<AppProps> = () => {
  const [vusList, setVusList] = useState<IVus[]>(undefined);

  return (
    <div className={styles["container"]}>
      {/* <PublicationSearch publicationService={publicationService} /> */}
      {/* <VusFileUpload vusService={vusService} /> */}
      {vusList && <ViewAllVus vusList={vusList} />}
      <button onClick={handleVusLoad}>Click to load VUS</button>
    </div>
  );

  function handleVusLoad(e: any) {
    e.preventDefault();

    vusService?.loadAllVus().then((res) => {
      if (res.isSuccess) {
        setVusList(res.vusList);
      } else {
        //TODO: Handle error
      }
    });
  }
};

export default App;
