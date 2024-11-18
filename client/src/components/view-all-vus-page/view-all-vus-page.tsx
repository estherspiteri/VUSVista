import React, { useEffect, useState } from "react";
import styles from "./view-all-vus-page.module.scss";
import { VusService } from "../../services/vus/vus.service";
import Loader from "../../atoms/loader/loader";
import VusTable from "./vus-table/vus-table";
import { IVUSSummary } from "../../models/vus-summary.model";

type ViewAllVusProps = { vusService: VusService };

const ViewAllVus: React.FunctionComponent<ViewAllVusProps> = (
  props: ViewAllVusProps
) => {
  const [vusList, setVusList] = useState<IVUSSummary[]>(undefined);

  useEffect(() => {
    props.vusService?.loadAllVus().then((res) => {
      if (res.isSuccess) {
        setVusList(res.vusList);
      }
    });
  }, [props.vusService]);

  return (
    <div className={styles["view-all-vus-container"]}>
      <div className={styles.title}>VUS List</div>
      <div className={styles.description}>
        <p>
          Below you can find a list of all the VUS stored within our database.
          Multiple column sorting can be enabled by holding down the SHIFT
          button and clicking on the desired columns.
        </p>
      </div>
      {vusList ? <VusTable vusList={vusList} /> : <Loader />}
    </div>
  );
};

export default ViewAllVus;
