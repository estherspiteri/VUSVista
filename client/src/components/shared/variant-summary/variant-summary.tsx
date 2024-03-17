import React from "react";
import styles from "./variant-summary.module.scss";
import { IVUSSummary } from "../../../models/publication-view.model";

type VariantSummaryProps = {
  variant: IVUSSummary;
};

const VariantSummary: React.FunctionComponent<VariantSummaryProps> = (
  props: VariantSummaryProps
) => {
  return (
    <div className={styles["variant-summary-container"]}>
      <div>
        <p className={styles["detail-name"]}>Id</p>
        <p className={styles.detail}>{props.variant.id}</p>
      </div>
      <div>
        <p className={styles["detail-name"]}>Chromosome</p>
        <p className={styles.detail}>{props.variant.chromosome}</p>
      </div>
      <div>
        <p className={styles["detail-name"]}>Position</p>
        <p className={styles.detail}>{props.variant.chromosomePosition}</p>
      </div>
      <div>
        <p className={styles["detail-name"]}>Gene</p>
        <p className={styles.detail}>{props.variant.gene}</p>
      </div>
      <div>
        <p className={styles["detail-name"]}>Reference</p>
        <p className={styles.detail}>{props.variant.refAllele}</p>
      </div>
      <div>
        <p className={styles["detail-name"]}>Alternate</p>
        <p className={styles.detail}>{props.variant.altAllele}</p>
      </div>
    </div>
  );
};

export default VariantSummary;
