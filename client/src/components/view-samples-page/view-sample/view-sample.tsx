import React, { useEffect, useRef, useState } from "react";
import styles from "./view-sample.module.scss";
import { IVus } from "../../../models/view-vus.model.tsx/view-vus.model";

type ViewSampleProps = {
  sample: {
    id: string;
    dateCollected: string;
    description: string;
    genomeVersion: string;
    // variants: IVus[];
    fileName: string;
    numOfVariants: number;
  };
  isColoured: boolean;
  onClickCallback?: (sampleId: string) => void;
};

const ViewSample: React.FunctionComponent<ViewSampleProps> = (
  props: ViewSampleProps
) => {
  return (
    <div
      className={`${styles["view-vus-container"]} ${
        props.isColoured ? styles.coloured : ""
      }`}
    >
      <div
        className={styles.header}
        onClick={() =>
          props.onClickCallback && props.onClickCallback(props.sample.id)
        }
      >
        <div className={styles["header-content"]}>{props.sample.id}</div>
        <div className={styles["header-content"]}>
          {props.sample.numOfVariants}
        </div>
        <div className={styles["header-content"]}>
          {props.sample.dateCollected}
        </div>
      </div>
    </div>
  );
};

export default ViewSample;
