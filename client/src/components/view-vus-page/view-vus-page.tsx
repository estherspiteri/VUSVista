import React, { useEffect, useState } from "react";
import styles from "./view-vus-page.module.scss";
import { IVus } from "../../models/view-vus.model.tsx/view-vus.model";
import ViewVus from "./view-vus/view-vus";
import { VusService } from "../../services/vus/vus.service";
import Loader from "../loader/loader";

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
      <div className={styles.title}>VUS LIST</div>
      <div className={styles.description}>
        <p>
          Below you can find a list of all the VUS stored within our database.
        </p>
      </div>
      {vusList ? (
        <div className={styles["view-all-vus"]}>
          <div className={styles.header}>
            <div>Chromosome</div>
            <div>Position</div>
            <div>Gene</div>
            <div>Reference</div>
            <div>Observed</div>
            <div>Genotype</div>
            <div>RSID</div>
            <div></div>
          </div>
          {vusList.map((vus, index) => {
            return <ViewVus vus={vus} isColoured={index % 2 === 0} />;
          })}
        </div>
      ) : (
        <Loader />
      )}
    </div>
  );
};

export default ViewAllVus;
