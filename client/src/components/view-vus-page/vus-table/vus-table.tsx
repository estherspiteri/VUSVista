import React from "react";
import styles from "./vus-table.module.scss";
import { IVus } from "../../../models/view-vus.model";
import ViewVus from "../../view-vus-page/view-vus/view-vus";

type VusTableProps = { vusList: IVus[]; showGenotype?: boolean };

const VusTable: React.FunctionComponent<VusTableProps> = (
  props: VusTableProps
) => {
  return (
    <div className={styles["vus-table-container"]}>
      <div
        className={`${styles.header} ${
          props.showGenotype ? styles["genotype-included"] : ""
        }`}
      >
        <div>Chromosome</div>
        <div>Position</div>
        <div>Gene</div>
        <div>Reference</div>
        <div>Observed</div>
        {props.showGenotype && <div>Genotype</div>}
        <div>RSID</div>
        <div></div>
      </div>
      {props.vusList.map((vus, index) => {
        return (
          <ViewVus
            vus={vus}
            isColoured={index % 2 === 0}
            showGenotype={props.showGenotype}
          />
        );
      })}
    </div>
  );
};

VusTable.defaultProps = {
  showGenotype: false,
};

export default VusTable;
