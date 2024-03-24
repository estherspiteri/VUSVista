import React, { useEffect, useRef, useState } from "react";
import styles from "./view-sample.module.scss";
import { IVus } from "../../../models/view-vus.model";
import { ISample, ISampleSummary } from "../../../models/view-samples.model";
import { Link } from "react-router-dom";

type ViewSampleProps = {
  sample: ISampleSummary;
  isColoured: boolean;
  onClickCallback?: (sampleId: string) => void;
};

const ViewSample: React.FunctionComponent<ViewSampleProps> = (
  props: ViewSampleProps
) => {
  return (
    <Link
      to={`/sample/${props.sample.sampleId}`}
      className={`${styles["view-vus-container"]} ${
        props.isColoured ? styles.coloured : ""
      } `}
    >
      <div
        className={styles.header}
        onClick={() =>
          props.onClickCallback && props.onClickCallback(props.sample.sampleId)
        }
      >
        <div className={`${styles["header-content"]} ${styles.id}`}>
          {props.sample.sampleId}
        </div>
        <div className={`${styles["header-content"]} ${styles.variants}`}>
          {props.sample.numOfVariants}
        </div>
      </div>
    </Link>
  );
};

export default ViewSample;
