import React, { useState } from "react";
import styles from "./view-vus-page.module.scss";
import { IVus } from "../../models/view-vus.model.tsx/view-vus.model";
import ViewAllVus from "./view-all-vus/view-all-vus";
import { VusService } from "../../services/vus/vus.service";

type ViewVusPageProps = {
  vusService?: VusService;
};

const ViewVusPage: React.FunctionComponent<ViewVusPageProps> = (
  props: ViewVusPageProps
) => {
  const [vusList, setVusList] = useState<IVus[]>(undefined);

  return (
    <div className={styles["view-vus-page-container"]}>
      {vusList && <ViewAllVus vusList={vusList} />}
      <button onClick={handleVusLoad}>Click to load VUS</button>
    </div>
  );

  function handleVusLoad(e: any) {
    e.preventDefault();

    props.vusService?.loadAllVus().then((res) => {
      if (res.isSuccess) {
        setVusList(res.vusList);
      } else {
        //TODO: Handle error
      }
    });
  }
};

export default ViewVusPage;
