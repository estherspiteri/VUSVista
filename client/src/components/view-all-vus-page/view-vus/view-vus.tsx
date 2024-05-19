import React from "react";
import styles from "./view-vus.module.scss";
import { Link } from "react-router-dom";
import { IVUSSummary } from "../../../models/vus-summary.model";
import { Row, flexRender } from "@tanstack/react-table";

type ViewVusProps = {
  vusRow?: Row<IVUSSummary>;
  isColoured?: boolean;
};

const ViewVus: React.FunctionComponent<ViewVusProps> = (
  props: ViewVusProps
) => {
  return (
    <Link
      to={`/vus/${props.vusRow.id}`}
      className={`${styles["view-vus-container"]} ${
        props.isColoured ? styles.coloured : ""
      }`}
    >
      <tr key={props.vusRow.id} className={styles.header}>
        {props.vusRow.getVisibleCells().map((cell) => (
          <td key={cell.id} className={styles["header-content"]}>
            {flexRender(cell.column.columnDef.cell, cell.getContext())}
          </td>
        ))}
      </tr>
      <div className={styles["additional-info"]}>
        <div className={styles["additional-info-content"]}></div>
      </div>
    </Link>
  );
};

export default ViewVus;
