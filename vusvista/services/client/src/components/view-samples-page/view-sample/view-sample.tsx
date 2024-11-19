import React from "react";
import styles from "./view-sample.module.scss";
import { ISampleSummary } from "../../../models/view-samples.model";
import { Link } from "react-router-dom";
import { Row, flexRender } from "@tanstack/react-table";

type ViewSampleProps = {
  sampleRow: Row<ISampleSummary>;
  isColoured: boolean;
  isClickable?: boolean;
  showCheckbox?: boolean;
  onCheckboxToggle?: () => void;
};

const ViewSample: React.FunctionComponent<ViewSampleProps> = (
  props: ViewSampleProps
) => {
  const row = (
    <tr key={props.sampleRow.id} className={styles.header}>
      {props.sampleRow.getVisibleCells().map((cell) => (
        <td
          key={cell.id}
          className={`${styles["header-content"]} ${
            cell.column.id === "sampleId" ? styles.id : styles.variants
          }`}
        >
          {flexRender(cell.column.columnDef.cell, cell.getContext())}
        </td>
      ))}
      {props.showCheckbox && (
        <input
          type="checkbox"
          className={styles.checkbox}
          name={props.sampleRow.original.sampleId}
          onChange={() => props.onCheckboxToggle && props.onCheckboxToggle()}
        />
      )}
    </tr>
  );

  const className = `${styles["view-sample-container"]} ${
    props.isColoured ? styles.coloured : ""
  }`;

  if (props.isClickable) {
    return (
      <Link
        to={`/sample/${props.sampleRow.original.sampleId}`}
        className={`${className} ${styles.clickable}`}      >
        {row}
      </Link>
    );
  } else {
    return <div className={className}>{row}</div>;
  }
};

ViewSample.defaultProps = {
  isClickable: true,
  showCheckbox: false,
};

export default ViewSample;
