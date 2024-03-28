import React from "react";
import styles from "./view-vus.module.scss";
import { Link } from "react-router-dom";
import { IVUSSummary } from "../../../models/vus-summary.model";

type ViewVusProps = {
  vus?: IVUSSummary;
  isColoured?: boolean;
};

const ViewVus: React.FunctionComponent<ViewVusProps> = (
  props: ViewVusProps
) => {
  return (
    <Link
      to={`/vus/${props.vus.id}`}
      className={`${styles["view-vus-container"]} ${
        props.isColoured ? styles.coloured : ""
      }`}
    >
      <div className={styles.header}>
        <div className={styles["header-content"]}>{props.vus.id}</div>
        <div className={styles["header-content"]}>{props.vus.chromosome}</div>
        <div className={styles["header-content"]}>
          {props.vus.chromosomePosition}
        </div>
        <div className={styles["header-content"]}>{props.vus.gene}</div>
        <div className={styles["header-content"]}>{props.vus.refAllele}</div>
        <div className={styles["header-content"]}>{props.vus.altAllele}</div>
        <div className={styles["header-content"]}>
          {props.vus.rsidDbsnpVerified ? props.vus.rsid : "-"}
        </div>
      </div>
      <div className={styles["additional-info"]}>
        <div className={styles["additional-info-content"]}></div>
      </div>
    </Link>
  );
};

export default ViewVus;
