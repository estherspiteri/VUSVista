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
            className={`${styles.pill} ${styles.clinvar} ${
              props.vus.clinvarErrorMsg.length > 0 ? styles.disabled : ""
            }`}
          >
            ClinVar
          </div>
          <div
            className={`${styles.pill} ${styles.dbsnp} ${
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
        {/*TODO: Next to is RSID verified do info icon - on hover show what info was compared. Same for clinvar*/}
        <div className={styles["additional-info-content"]}>
          <div
            className={`${styles["info-container"]} ${styles["clinvar-info"]} ${
              props.vus.clinvarErrorMsg.length > 0 ? styles.disabled : ""
            }`}
          >
            <p className={styles["info-cotainer-title"]}>ClinVar</p>
            <div className={styles.info}>
              {props.vus.clinvarClassification.length > 0 ? (
                <>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Classification:</div>
                    {props.vus.clinvarClassification}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Review status:</div>
                    {props.vus.clinvarClassificationReviewStatus}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Last evaluated:</div>
                    {props.vus.clinvarClassificationLastEval}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Canonical SPDI:</div>
                    {props.vus.clinvarCanonicalSpdi}
                  </div>
                </>
              ) : (
                <div className={styles.information}>
                  <div className={styles["info-title"]}>Error message:</div>
                  {props.vus.clinvarErrorMsg}
                </div>
              )}
            </div>
          </div>
          <div
            className={`${styles["info-container"]} ${styles["dbsnp-info"]} ${
              props.vus.rsidDbsnpVerified ? "" : styles.disabled
            }`}
          >
            <p className={styles["info-container-title"]}>dbSnp</p>
            <div className={styles.info}>
              <div className={styles.information}>
                <div className={styles["info-title"]}>Is RSID verified:</div>
                {props.vus.rsidDbsnpVerified.toString()}
              </div>
              {!props.vus.rsidDbsnpVerified && (
                <div className={styles.information}>
                  <div className={styles["info-title"]}>Error message:</div>
                  {props.vus.rsidDbsnpErrorMsgs}
                </div>
              )}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ViewVus;
