import React from "react";
import styles from "./sample-info.module.scss";
import { Genotype, ISample } from "../../../models/view-samples.model";
import { SampleService } from "../../../services/sample/sample.service";
import VariantSummary from "../../shared/variant-summary/variant-summary";
import SamplePhenotypeSelection from "../sample-phenotype-selection/sample-phenotype-selection";
import Icon from "../../../atoms/icons/icon";
import { openInNewWindow } from "../../../helpers/open-links";

type SampleInfoProps = {
  sample: ISample;
  sampleService: SampleService;
};

const SampleInfo: React.FunctionComponent<SampleInfoProps> = (
  props: SampleInfoProps
) => {
  return (
    <div className={styles["sample-info-container"]}>
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
                <div className={styles["variant-summary"]}>Variant Summary</div>
                <div className={styles.genotype}>Genotype</div>
                <div>HGVS</div>
              </div>
              <div className={styles.variants}>
                {props.sample.variants.map((v) => {
                  return (
                    <div className={styles.variant}>
                      <Icon
                        name="external-link"
                        className={styles["external-link"]}
                        onClick={() => openInNewWindow(`/vus/${v.variantId}`)}
                      />
                      <div className={styles["variant-info-container"]}>
                        <div className={styles["variant-info"]}>
                          <div className={styles["variant-summary"]}>
                            {/*TODO: on click take to VUS page*/}
                            <VariantSummary variant={v.variant} />
                          </div>
                          <div className={styles.genotype}>
                            {v.genotype === Genotype.Heterozygous ? "Aa" : "AA"}
                          </div>
                          <div>{v.hgvs}</div>
                        </div>
                      </div>
                    </div>
                  );
                })}
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default SampleInfo;
