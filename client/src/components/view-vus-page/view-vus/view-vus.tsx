import React, { useEffect, useRef, useState } from "react";
import styles from "./view-vus.module.scss";
import { IVus } from "../../../models/view-vus.model.tsx/view-vus.model";

type ViewVusProps = {
  vus?: IVus;
  isColoured?: boolean;
};

const ViewVus: React.FunctionComponent<ViewVusProps> = (
  props: ViewVusProps
) => {
  const ref = useRef<HTMLDivElement>(null);
  const [isAdditionalInfoVisible, setIsAdditionalInfoVisible] = useState(false);

  //close additional info on outside click
  useEffect(() => {
    function handleClickOutside(event) {
      if (ref.current && !ref.current.contains(event.target)) {
        setIsAdditionalInfoVisible(false);
      }
    }
    // Bind the event listener
    document.addEventListener("mousedown", handleClickOutside);

    return () => {
      // Unbind the event listener on clean up
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [ref]);

  return (
    <div
      ref={ref}
      className={`${styles["view-vus-container"]} ${
        props.isColoured ? styles.coloured : ""
      }`}
    >
      <div
        className={`${styles.header}  `}
        onClick={() => setIsAdditionalInfoVisible(!isAdditionalInfoVisible)}
      >
        <div className={styles["header-content"]}>{props.vus.chromosome}</div>
        <div className={styles["header-content"]}>
          {props.vus.chromosomePosition}
        </div>
        <div className={styles["header-content"]}>{props.vus.gene}</div>
        <div className={styles["header-content"]}>{props.vus.refAllele}</div>
        <div className={styles["header-content"]}>
          {props.vus.observedAllele}
        </div>
        <div className={styles["header-content"]}>{props.vus.genotype}</div>
        <div className={styles["header-content"]}>{props.vus.rsid}</div>
        <div className={styles.pills}>
          <div
            className={`${styles.clinvar} ${
              props.vus.clinvarErrorMsg.length > 0 ? styles.disabled : ""
            }`}
          >
            ClinVar
          </div>
          <div
            className={`${styles.dbsnp} ${
              props.vus.rsidDbsnpVerified ? "" : styles.disabled
            }`}
          >
            dbSNP
          </div>
        </div>
      </div>
      <div
        className={`${styles["additional-info"]} ${
          isAdditionalInfoVisible ? styles["visible-additional-info"] : ""
        }`}
      >
        <div className={styles["additional-info-content"]}>
          {props.vus.clinvarClassification}
          &nbsp;
          {props.vus.clinvarErrorMsg}
        </div>
      </div>
    </div>
  );
};

export default ViewVus;
