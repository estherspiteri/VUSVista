import React from "react";
import styles from "./view-sample.module.scss";
import { ISampleSummary } from "../../../models/view-samples.model";
import { Link } from "react-router-dom";
import { Row, flexRender } from "@tanstack/react-table";

type ViewSampleProps = {
  sampleRow: Row<ISampleSummary>;
  // sample: ISampleSummary;
  isColoured: boolean;
};

const ViewSample: React.FunctionComponent<ViewSampleProps> = (
  props: ViewSampleProps
) => {
  return (
    <Link
      to={`/sample/${props.sampleRow.original.sampleId}`}
      className={`${styles["view-vus-container"]} ${
        props.isColoured ? styles.coloured : ""
      } `}
    >
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
      </tr>
    </Link>
  );
};

export default ViewSample;
