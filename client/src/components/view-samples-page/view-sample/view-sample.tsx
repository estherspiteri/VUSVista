import React, { useEffect, useRef, useState } from "react";
import styles from "./view-sample.module.scss";
import { IVus } from "../../../models/view-vus.model";
import { ISample } from "../../../models/view-samples.model";

type ViewSampleProps = {
  sample: ISample;
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
          props.onClickCallback && props.onClickCallback(props.sample.sampleId)
        }
      >
        <div className={`${styles["header-content"]} ${styles.id}`}>{props.sample.sampleId}</div>
        <div className={`${styles["header-content"]} ${styles.variants}`}>
          {props.sample.variants.length}
        </div>
        <div className={`${styles["header-content"]} ${styles.date}`}>
          {`${props.sample.dateOfFileUpload.getDate()}/${
            props.sample.dateOfFileUpload.getMonth() + 1
          }/${props.sample.dateOfFileUpload.getFullYear()}`}
        </div>
      </div>
    </div>
  );
};

export default ViewSample;
