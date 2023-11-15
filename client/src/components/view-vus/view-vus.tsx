import React, { useRef, useState } from "react";
import styles from "./view-vus.module.scss";
import { IVus } from "../../models/view-vus.model.tsx/view-vus.model";

type ViewVusProps = {
  vusList?: IVus[];
};

const ViewVus: React.FunctionComponent<ViewVusProps> = (
  props: ViewVusProps
) => {
  return (
    <div className={styles["view-vus-container"]}>
      {props.vusList.map((vus) => (
        <p>
          {vus.vusID} {vus.rsid} {vus.chromosome} {vus.chromosomePosition}{" "}
          {vus.gene} {vus.refAllele} {vus.observedAllele} {vus.genotype}{" "}
          {vus.rsid} {vus.rsidDbsnpVerified} {vus.clinvarClassification}{" "}
          {vus.clinvarErrorMsg}
        </p>
      ))}
    </div>
  );
};

export default ViewVus;
