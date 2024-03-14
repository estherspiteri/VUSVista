import React from "react";
import styles from "./sample-info.module.scss";
import { Genotype, ISample } from "../../../models/view-samples.model";
import Icon from "../../../atoms/icon/icon";
import AcmgRulesEdit from "../acmg-rules-edit/acmg-rules-edit";
import { SampleService } from "../../../services/sample/sample.service";

type SampleInfoProps = {
  sample: ISample;
  acmgRuleNames: string[];
  sampleService: SampleService
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
            <div className={`${styles.info} ${styles["variant-phenotypes"]}`}>
              {props.sample.phenotype.map((p) => {
                return (
                  <div className={styles.phenotype}>
                    <div>&#x2022;</div>&nbsp;
                    <div>
                      <span className={styles["ontology-id"]}>
                        {p.ontologyId}
                      </span>
                      : {p.name}
                    </div>
                  </div>
                );
              })}
              {/* </div> */}
            </div>
          </div>

          <div className={`${styles.info} ${styles["info-container"]}`}>
            <div className={styles["info-title"]}>Variants:</div>
            <div className={styles["info-container-content"]}>
              <div className={styles["variant-titles"]}>
                <div className={styles["variant-id"]}>Variant Id</div>
                <div className={styles.genotype}>Genotype</div>
                <div className={styles.acmg}>Assigned ACMG rules</div>
              </div>
              <div className={styles.variants}>
                {props.sample.variants.map((v) => {
                  return (
                    <div className={styles.variant}>
                      <div className={styles["variant-id"]}>{v.variantId}</div>
                      <div className={styles.genotype}>
                        {v.genotype == Genotype.Heterozygous ? "Aa" : "AA"}
                      </div>
                      <AcmgRulesEdit
                        variantId={v.variantId}
                        sampleId={props.sample.sampleId}
                        variantAcmgRules={v.acmgRuleNames}
                        allAcmgRules={props.acmgRuleNames}
                        sampleService={props.sampleService}
                      />
                    </div>
                  );
                })}
              </div>
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
