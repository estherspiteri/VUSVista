import React from "react";
import styles from "./sample-info.module.scss";
import { ISample } from "../../../models/view-samples.model";
import Icon from "../../../atoms/icon/icon";

type SampleInfoProps = {
  sample: ISample;
};

const SampleInfo: React.FunctionComponent<SampleInfoProps> = (
  props: SampleInfoProps
) => {
  return (
    <div
      className={`${styles["sample-info-container"]} ${
        props.sample ? "" : styles["no-sample-selected-container"]
      }`}
    >
      <p className={styles["sample-info-title"]}>Sample information</p>
      {props.sample ? (
        <div className={styles["sample-info"]}>
          <div className={styles.info}>
            <div className={styles["info-title"]}>Sample Id:</div>
            {props.sample.sampleId}
          </div>

          <div className={styles.info}>
            <div className={styles["info-title"]}>Genome Version:</div>
            {props.sample.genomeVersion}
          </div>

          <div className={styles.info}>
            <div className={styles["info-title"]}>File upload name:</div>
            {props.sample.fileUploadName}
          </div>

          <div className={`${styles.info} ${styles["info-container"]}`}>
            <div className={styles["info-title"]}>Phenotypes:</div>
            <div className={styles["info-container-content"]}>
              <div className={`${styles.info} ${styles["variant-phenotypes"]}`}>
                {props.sample.phenotype.map((p) => {
                  return (
                    <div className={styles.phenotype}>
                      <div>&#x2022;</div>&nbsp;
                      <div>
                        <b>{p.ontologyId}</b>: {p.name}
                      </div>
                    </div>
                  );
                })}
              </div>
            </div>
          </div>

          <div className={`${styles.info} ${styles["info-container"]}`}>
            <div className={styles["info-title"]}>Variants:</div>
            <div className={styles["info-container-content"]}>
              <div className={styles["variant-titles"]}>
                <div className={styles["variant-id"]}>Variant Ids</div>
                <div>Genotypes</div>
                {/* <div>Assigned ACMG rules</div> */}
              </div>
              {props.sample.variants.map((v) => {
                return (
                  <div className={styles.variant}>
                    <div className={styles["variant-id"]}>{v.variantId}</div>
                    <div>{v.genotype}</div>
                    {/* <div className={styles.rules}>
                <div className={styles["pathogenic-supporting"]}>
                  PP1
                </div>
                <div className={styles["pathogenic-medium"]}>PM6</div>
                <div className={styles["pathogenic-strong"]}>PS2</div>
              </div> */}
                  </div>
                );
              })}
            </div>
          </div>
        </div>
      ) : (
        <div className={styles["no-selected-sample"]}>
          <Icon name="left-arrow-circle" stroke="#fff" />
          <span>
            Select a Sample from the left to view information related to it.
          </span>
        </div>
      )}
    </div>
  );
};

export default SampleInfo;
