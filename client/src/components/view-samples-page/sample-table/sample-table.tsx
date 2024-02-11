import React from "react";
import styles from "./sample-table.module.scss";
import ViewSample from "../view-sample/view-sample";
import { ISample } from "../../../models/view-samples.model";

type SampleTableProps = {
  sampleList: ISample[];
  onSampleClickCallback?: (sampleId: string) => void;
};

const SampleTable: React.FunctionComponent<SampleTableProps> = (
  props: SampleTableProps
) => {
  return (
    <div className={styles["sample-table-container"]}>
      <div className={styles.header}>
        <div className={`${styles['header-content']} ${styles.id}`}>Sample Id</div>
        <div className={`${styles['header-content']} ${styles.variants}`}>No. of Variants</div>
        <div className={`${styles['header-content']} ${styles.date}`}>Date uploaded</div>
      </div>
      {props.sampleList.map((sample, index) => {
        return (
          <ViewSample
            sample={sample}
            isColoured={index % 2 === 0}
            onClickCallback={props.onSampleClickCallback}
          />
        );
      })}
    </div>
  );
};

export default SampleTable;
