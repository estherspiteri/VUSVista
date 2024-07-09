import React from "react";
import styles from "./view-vus.module.scss";
import { Link } from "react-router-dom";
import { IVUSSummary } from "../../../models/vus-summary.model";
import { Row, flexRender } from "@tanstack/react-table";

type ViewVusProps = {
  vusRow?: Row<IVUSSummary>;
  isColoured?: boolean;
  isClickable?: boolean;
  showCheckbox?: boolean;
  onCheckboxToggle?: () => void;
  onVariantClickCallback?: () => void;
};

const ViewVus: React.FunctionComponent<ViewVusProps> = (
  props: ViewVusProps
) => {
  const row = (
    <tr key={props.vusRow.id} className={styles.header}>
      {props.vusRow.getVisibleCells().map((cell) => (
        <td key={cell.id} className={styles["header-content"]}>
          {flexRender(cell.column.columnDef.cell, cell.getContext())}
        </td>
      ))}
      {props.showCheckbox && (
        <input
          type="checkbox"
          className={styles.checkbox}
          name={props.vusRow.original.id.toString()}
          onChange={() => props.onCheckboxToggle && props.onCheckboxToggle()}
        />
      )}
    </tr>
  );

  const className = `${styles["view-vus-container"]} ${
    props.isColoured ? styles.coloured : ""
  }`;

  if (props.isClickable) {
    return (
      <Link
        to={`/vus/${props.vusRow.original.id}`}
        onClick={() =>
          props.onVariantClickCallback && props.onVariantClickCallback()
        }
        className={`${className} ${styles.clickable}`}
      >
        {row}
      </Link>
    );
  } else {
    return <div className={className}>{row}</div>;
  }
};

ViewVus.defaultProps = {
  isClickable: true,
  showCheckbox: false,
};

export default ViewVus;
