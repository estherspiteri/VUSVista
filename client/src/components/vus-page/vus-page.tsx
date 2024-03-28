import React from "react";
import styles from "./vus-page.module.scss";
import { IVus } from "../../models/view-vus.model";
import VusInfo from "./vus-info/vus-info";

type VusPageProps = {
  vus: IVus;
};

const VusPage: React.FunctionComponent<VusPageProps> = (
  props: VusPageProps
) => {
  return (
    <div className={styles["sample-page-container"]}>
      <div className={styles.title}>Vus Information</div>
      <div className={styles.description}>
        <p>Below you can find information about the selected variant.</p>
      </div>

      <VusInfo
        vus={props.vus}
      />
    </div>
  );
};

export default VusPage;
