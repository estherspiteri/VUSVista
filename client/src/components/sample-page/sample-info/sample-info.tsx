import React, { useState } from "react";
import styles from "./sample-info.module.scss";
import { Genotype, ISample } from "../../../models/view-samples.model";
import Icon from "../../../atoms/icon/icon";
import AcmgRulesEdit from "../acmg-rules-edit/acmg-rules-edit";
import { SampleService } from "../../../services/sample/sample.service";
import { IAcmgRule } from "../../../models/acmg-rule.model";
import VariantSummary from "../../shared/variant-summary/variant-summary";
import AcmgRuleInfo from "../acmg-rule-info/acmg-rule-info";
import SamplePhenotypeSelection from "../sample-phenotype-selection/sample-phenotype-selection";

type SampleInfoProps = {
  sample: ISample;
  acmgRules: IAcmgRule[];
  sampleService: SampleService;
};

const SampleInfo: React.FunctionComponent<SampleInfoProps> = (
  props: SampleInfoProps
) => {
  const [acmgRuleHover, setAcmgRuleHover] = useState<number | undefined>(
    undefined
  );

  return (
    <div
      className={`${styles["sample-info-container"]} ${
        props.sample ? "" : styles["no-sample-selected-container"]
      }`}
    >
      {props.sample ? (
        <div className={styles["sample-info"]}>
          <div className={styles["top-container"]}>
            {/** General information */}
            <div className={styles.information}>
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
            </div>
          </div>

          {/** Phenotypes */}
          <div className={styles.phenotypes}>
            <div className={styles["info-title"]}>Phenotypes:</div>
            <SamplePhenotypeSelection
              sampleService={props.sampleService}
              sampleId={props.sample.sampleId}
              selectedPhenotypes={props.sample.phenotype}
            />
          </div>

          {/** Variants */}
          <div className={styles["sample-variant-container"]}>
            <div className={styles["info-title"]}>Variants:</div>
            <div className={styles["sample-variant-container-content"]}>
              <div className={styles["variants-container"]}>
                <div className={styles["variant-titles"]}>
                  <div className={styles["variant-summary"]}>
                    Variant Summary
                  </div>
                  <div className={styles.genotype}>Genotype</div>
                  <div className={styles.acmg}>Assigned ACMG rules</div>
                </div>
                <div className={styles.variants}>
                  {props.sample.variants.map((v) => {
                    return (
                      <div className={styles.variant}>
                        <div className={styles["variant-summary"]}>
                          {/*TODO: on click take to VUS page*/}
                          <VariantSummary variant={v.variant} />
                        </div>
                        <div className={styles.genotype}>
                          {v.genotype === Genotype.Heterozygous ? "Aa" : "AA"}
                        </div>
                        <AcmgRulesEdit
                          variantId={v.variantId}
                          sampleId={props.sample.sampleId}
                          variantAcmgRuleIds={v.acmgRuleIds}
                          allAcmgRules={props.acmgRules}
                          sampleService={props.sampleService}
                          onMenuAcmgRuleHover={(acmgRuleId?: number) =>
                            setAcmgRuleHover(acmgRuleId)
                          }
                        />
                      </div>
                    );
                  })}
                </div>
              </div>
              <div className={styles["acmg-info"]}>
                <AcmgRuleInfo
                  acmgRule={props.acmgRules.find((r) => r.id === acmgRuleHover)}
                />
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
