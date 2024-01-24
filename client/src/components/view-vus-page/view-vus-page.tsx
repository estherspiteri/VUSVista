import React, { useEffect, useState } from "react";
import styles from "./view-vus-page.module.scss";
import { IVus } from "../../models/view-vus.model";
import { VusService } from "../../services/vus/vus.service";
import Loader from "../loader/loader";
import VusTable from "./vus-table/vus-table";

type ViewAllVusProps = { vusService: VusService };

const ViewAllVus: React.FunctionComponent<ViewAllVusProps> = (
  props: ViewAllVusProps
) => {
  const [vusList, setVusList] = useState<IVus[]>(undefined);

  useEffect(() => {
    props.vusService?.loadAllVus().then((res) => {
      if (res.isSuccess) {
        setVusList(res.vusList);
      } else {
        //TODO: Handle error
      }
    });
  }, []);

  return (
    <div className={styles["view-all-vus-container"]}>
      <div className={styles.title}>VUS List</div>
      <div className={styles.description}>
        <p>
          Below you can find a list of all the VUS stored within our database.
        </p>
      </div>
      {vusList ? (
        <VusTable vusList={vusList} showGenotype={false} showZygosity={true} />
      ) : (
        <Loader />
      )}
    </div>
  );
};

export default ViewAllVus;
