import React from "react";
import styles from "./vus-table.module.scss";
import ViewVus from "../view-vus/view-vus";
import { IVUSSummary } from "../../../models/vus-summary.model";

type VusTableProps = {
  vusList: IVUSSummary[];
};

const VusTable: React.FunctionComponent<VusTableProps> = (
  props: VusTableProps
) => {
  return (
    <div className={styles["vus-table-container"]}>
      <div className={styles.header}>
        <div>Variant Id</div>
        <div>Chromosome</div>
        <div>Position</div>
        <div>Gene</div>
        <div>Reference</div>
        <div>Alternate</div>
        <div>RSID</div>
      </div>
      {props.vusList.map((vus, index) => {
        return <ViewVus vus={vus} isColoured={index % 2 === 0} />;
      })}
    </div>
  );
};

export default VusTable;
