import React from "react";
import styles from "./sample-table.module.scss";
import ViewSample from "../view-sample/view-sample";

type SampleTableProps = {
  sampleList: {
    id: string;
    dateCollected: string;
    description: string;
    genomeVersion: string;
    // variants: IVus[];
    fileName: string;
    numOfVariants: number;
  }[];
  onSampleClickCallback?: (sampleId: string) => void;
};

const SampleTable: React.FunctionComponent<SampleTableProps> = (
  props: SampleTableProps
) => {
  return (
    <div className={styles["sample-table-container"]}>
      <div className={styles.header}>
        <div>Sample Id</div>
        <div>No. of Variants</div>
        <div>Date collected</div>
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
