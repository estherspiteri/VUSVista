import React, { useEffect, useRef, useState } from "react";
import styles from "./view-vus.module.scss";
import { IVus } from "../../../models/view-vus.model";
import Icon from "../../../atoms/icons/icon";
import { Link } from "react-router-dom";

type ViewVusProps = {
  vus?: IVus;
  isColoured?: boolean;
  showGenotype?: boolean;
  showZygosityQty?: boolean;
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
      <div
        className={`${styles.header} ${
          props.showGenotype ? styles["genotype-included"] : ""
        }`}
      >
        <div className={styles["header-content"]}>{props.vus.id}</div>
        <div className={styles["header-content"]}>{props.vus.chromosome}</div>
        <div className={styles["header-content"]}>
          {props.vus.chromosomePosition}
        </div>
        <div className={styles["header-content"]}>{props.vus.gene}</div>
        <div className={styles["header-content"]}>{props.vus.refAllele}</div>
        <div className={styles["header-content"]}>{props.vus.altAllele}</div>
        {props.showGenotype && (
          <div className={styles["header-content"]}>{props.vus.genotype}</div>
        )}
        <div className={styles["header-content"]}>
          {props.vus.rsidDbsnpVerified ? props.vus.rsid : "-"}
        </div>
        {/* <div className={styles.pills}>
          <div
            className={`${styles.pill} ${styles.dbsnp} ${
              props.vus.rsidDbsnpVerified
                ? styles.active
                : props.vus.rsid.length > 0
                ? styles["unverified-rsid"]
                : styles.disabled
            }`}
            onClick={(e) => {
              if (props.vus.rsidDbsnpVerified) {
                e.stopPropagation();
                openInNewWindow(
                  `https://www.ncbi.nlm.nih.gov/snp/${props.vus.rsid}`
                );
              }
            }}
          >
            dbSNP
          </div>
          <div
            className={`${styles.pill} ${styles.clinvar} ${
              props.vus.clinvarClassification.length === 0
                ? styles.disabled
                : props.vus.rsid.length > 0 && !props.vus.rsidDbsnpVerified
                ? styles["unverified-rsid"]
                : styles.active
            }`}
            onClick={(e) => {
              if (props.vus.clinvarErrorMsg.length === 0) {
                e.stopPropagation();
                openInNewWindow(
                  `https://www.ncbi.nlm.nih.gov/clinvar/variation/${props.vus.clinvarUid}`
                );
              }
            }}
          >
            ClinVar
          </div>
          {props.vus.rsidDbsnpVerified && (
            <Link
              to={`/publication-view/${props.vus.variantId}`}
              className={styles["pub-icon"]}
            >
              <Icon name="publication" />
            </Link>
          )}
        </div> */}
      </div>
      <div className={styles["additional-info"]}>
        {/*TODO: Next to is RSID verified do info icon - on hover show what info was compared. Same for clinvar*/}
        {/* <div className={styles["additional-info-content"]}>
          {props.showZygosityQty && (
            <>
              <p>Homozygotes: {props.vus.numHomozygous}</p>
              <p>Heterozygotes: {props.vus.numHeterozygous}</p>
            </>
          )}
          <div
            className={`${styles["info-container"]} ${styles["dbsnp-info"]} ${
              props.vus.rsidDbsnpVerified
                ? ""
                : props.vus.rsid.length > 0
                ? styles["unverified-rsid"]
                : styles.disabled
            }`}
          >
            <p className={styles["info-container-title"]}>dbSnp</p>

            <div className={styles.info}>
              {props.vus.rsid.length > 0 ? (
                <>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>
                      Is RSID verified:
                    </div>
                    {props.vus.rsidDbsnpVerified.toString()}
                  </div>
                  {!props.vus.rsidDbsnpVerified && (
                    <>
                      {props.vus.rsid !== "NORSID" && (
                        <div className={styles.information}>
                          <div className={styles["info-title"]}>
                            Suggested RSID:
                          </div>
                          <a
                            href={`https://www.ncbi.nlm.nih.gov/snp/${props.vus.rsid}`}
                            target="_blank"
                          >
                            {props.vus.rsid}
                          </a>
                        </div>
                      )}
                      <div className={styles.information}>
                        <div className={styles["info-title"]}>
                          Error message:
                        </div>
                        {props.vus.rsid === "NORSID"
                          ? "No RSID found."
                          : props.vus.rsidDbsnpErrorMsgs}
                      </div>
                    </>
                  )}
                </>
              ) : (
                <div className={styles.information}>
                  No valid RSID found for this variant.
                </div>
              )}
            </div>
          </div>
          <div
            className={`${styles["info-container"]} ${styles["clinvar-info"]} ${
              props.vus.clinvarClassification.length === 0
                ? styles.disabled
                : props.vus.rsid.length > 0 && !props.vus.rsidDbsnpVerified
                ? styles["unverified-rsid"]
                : ""
            }`}
          >
            <p className={styles["info-cotainer-title"]}>
              ClinVar{" "}
              {props.vus.clinvarClassification.length > 0 &&
                !props.vus.rsidDbsnpVerified &&
                "of suggested dbSNP RSID"}
            </p>
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
              ) : props.vus.clinvarErrorMsg.length > 0 ? (
                <div className={styles.information}>
                  <div className={styles["info-title"]}>Error message:</div>
                  {props.vus.clinvarErrorMsg}
                </div>
              ) : (
                <div className={styles.information}>
                  No Clinvar entry found based on dbSNP's RSID
                </div>
              )}
            </div>
          </div>
        </div> */}
      </div>
    </Link>
  );

  function openInNewWindow(url: string) {
    const newWindow = window.open(url, "_blank", "noopener,noreferrer");
    if (newWindow) newWindow.opener = null;
  }
};

ViewVus.defaultProps = {
  showGenotype: false,
  showZygosityQty: false,
};

export default ViewVus;
