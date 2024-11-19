import React, { useEffect, useState } from "react";
import { Column, RowData } from "@tanstack/react-table";
import DebouncedInput from "../debounced-input/debounced-input";
import styles from "./filter.module.scss";
import Icon from "../../../../atoms/icons/icon";

type FilterProps = {
  id: string;
  column: Column<any, unknown>;
};

declare module "@tanstack/react-table" {
  //allows us to define custom properties for our columns
  interface ColumnMeta<TData extends RowData, TValue> {
    filterVariant?: "text" | "number" | "boolean";
  }
}

const Filter: React.FunctionComponent<FilterProps> = (props: FilterProps) => {
  const columnFilterValue = props.column.getFilterValue();
  const { filterVariant } = props.column.columnDef.meta ?? {};

  const [checkboxState, setCheckboxState] = useState(0);

  useEffect(() => {
    if (checkboxState === 0) {
      props.column.setFilterValue(undefined);
    } else if (checkboxState === 1) {
      props.column.setFilterValue(true);
    } else {
      props.column.setFilterValue(false);
    }
  }, [checkboxState, props.column]);

  return filterVariant === "number" ? (
    <DebouncedInput
      onChange={(value) => props.column.setFilterValue(value)}
      placeholder={`Search...`}
      type="number"
      value={(columnFilterValue ?? "") as string}
      min={0}
    />
  ) : filterVariant === "boolean" ? (
    <div className={styles["checkbox-wrapper"]}>
      <div
        className={styles["checkbox-container"]}
        onClick={() =>
          setCheckboxState(checkboxState === 2 ? 0 : checkboxState + 1)
        }
      >
        {(checkboxState === 1 || checkboxState === 2) && (
          <Icon
            name={checkboxState === 1 ? "checkmark" : "close"}
            fill="#008080"
            stroke={checkboxState === 2 ? "#008080" : undefined}
          />
        )}
      </div>
    </div>
  ) : (
    <DebouncedInput
      onChange={(value) => props.column.setFilterValue(value)}
      placeholder={`Search...`}
      type="text"
      value={(columnFilterValue ?? "") as string}
    />
  );
};

export default Filter;
