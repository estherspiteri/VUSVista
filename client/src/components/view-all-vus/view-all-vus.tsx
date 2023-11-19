import React from "react";
import styles from "./view-all-vus.module.scss";
import { IVus } from "../../models/view-vus.model.tsx/view-vus.model";
import ViewVus from "../view-vus/view-vus";

type ViewAllVusProps = {
  vusList?: IVus[];
};

const ViewAllVus: React.FunctionComponent<ViewAllVusProps> = (
  props: ViewAllVusProps
) => {
  return (
    <div className={styles["view-all-vus-container"]}>
      <div className={styles.title}>VUS LIST</div>
      <div className={styles.description}>
        <p>
          Below you can find a list of all the VUS stored within our database.
        </p>
      </div>
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
        {props.vusList.map((vus, index) => {
          return <ViewVus vus={vus} isColoured={index % 2 === 0} />;
        })}
      </div>
    </div>
  );
};

export default ViewAllVus;
