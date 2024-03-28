import React from "react";
import styles from "./vus-info.module.scss";
import { IVus } from "../../../models/view-vus.model";
import VariantSummary from "../../shared/variant-summary/variant-summary";
import Icon from "../../../atoms/icons/icon";
import { openInNewWindow } from "../../../helpers/open-links";
import { Link } from "react-router-dom";

type VusInfoProps = {
  vus: IVus;
};

const VusInfo: React.FunctionComponent<VusInfoProps> = (
  props: VusInfoProps
) => {
  return (
    <div className={styles["vus-info-container"]}>
      <div className={styles["vus-info"]}>
        <div className={styles["vus-summary"]}>
          <VariantSummary
            variant={{
              id: props.vus.id,
              chromosome: props.vus.chromosome,
              chromosomePosition: props.vus.chromosomePosition,
              gene: props.vus.gene,
              refAllele: props.vus.refAllele,
              altAllele: props.vus.altAllele,
            }}
          />
        </div>
        <div className={styles["top-container"]}>
          {/** General information */}
          <div className={styles.information}>
            <div className={styles.info}>
              <div className={styles["info-title"]}>Type:</div>
              {props.vus.type}
            </div>

            <div className={styles.info}>
              <div className={styles["info-title"]}>Homozygotes:</div>
              {props.vus.numHomozygous}
            </div>

            <div className={styles.info}>
              <div className={styles["info-title"]}>Heterozygotes:</div>
              {props.vus.numHeterozygous}
            </div>

            <div className={styles.info}>
              <div className={styles["info-title"]}>Publications:</div>
              {props.vus.rsidDbsnpVerified && (
                <Icon
                  name="publication"
                  className={styles["pub-icon"]}
                  onClick={() =>
                    openInNewWindow(`/publication-view/${props.vus.id}`)
                  }
                />
              )}
            </div>
          </div>
        </div>

        {/*TODO: Next to is RSID verified do info icon - on hover show what info was compared. Same for clinvar*/}
        <div className={styles["external-ref"]}>
          {/** Clinvar */}
          <div className={styles["info-container"]}>
            <div className={styles["external-ref-title-container"]}>
              <p className={styles["info-title"]}>
                {`ClinVar${
                  props.vus.clinvarClassification?.length > 0 &&
                  !props.vus.rsidDbsnpVerified
                    ? " of suggested dbSNP RSID"
                    : ""
                }:`}
              </p>
              <Icon
                name="external-link"
                className={`${styles["external-link"]} ${styles.clinvar} ${
                  props.vus.clinvarClassification?.length === 0
                    ? styles.disabled
                    : props.vus.rsid.length > 0 && !props.vus.rsidDbsnpVerified
                    ? styles["unverified-rsid"]
                    : styles.active
                }`}
                onClick={(e) => {
                  if (props.vus.clinvarErrorMsg?.length === 0) {
                    e.stopPropagation();
                    openInNewWindow(
                      `https://www.ncbi.nlm.nih.gov/clinvar/variation/${props.vus.clinvarUid}`
                    );
                  }
                }}
              />
            </div>
            <div
              className={`${styles.info} ${
                props.vus.clinvarClassification?.length === 0
                  ? styles.disabled
                  : props.vus.rsid.length > 0 && !props.vus.rsidDbsnpVerified
                  ? styles["unverified-rsid"]
                  : ""
              }`}
            >
              {props.vus.clinvarClassification?.length > 0 ? (
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
              ) : props.vus.clinvarErrorMsg?.length > 0 ? (
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

          {/** DbSnp */}
          <div className={styles["info-container"]}>
            <div className={styles["external-ref-title-container"]}>
              <p className={styles["info-title"]}>dbSNP:</p>
              <Icon
                name="external-link"
                className={`${styles["external-link"]} ${styles.dbsnp} ${
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
              />
            </div>

            <div
              className={`${styles.info} ${
                props.vus.rsidDbsnpVerified
                  ? ""
                  : props.vus.rsid?.length > 0
                  ? styles["unverified-rsid"]
                  : styles.disabled
              }`}
            >
              {props.vus.rsid?.length > 0 ? (
                <>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>
                      Is RSID verified:
                    </div>
                    {props.vus.rsidDbsnpVerified.toString()}
                  </div>
                  {props.vus.rsidDbsnpVerified ? (
                    <div className={styles.information}>
                      <div className={styles["info-title"]}>RSID:</div>
                      {props.vus.rsid}
                    </div>
                  ) : (
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
        </div>

        <div className={styles["samples-container"]}>
          <p className={styles["info-title"]}>Samples with this variant:</p>
          <div className={styles.samples}>
            {props.vus.samples.map((s) => (
              <div className={styles.sample}>
                <Icon
                  name="external-link"
                  className={styles["external-link-sample"]}
                  onClick={() => openInNewWindow(`/sample/${s}`)}
                />
                <p>{s}</p>
              </div>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
};

export default VusInfo;
