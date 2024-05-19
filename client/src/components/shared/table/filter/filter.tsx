import React from "react";
import { Column, RowData } from "@tanstack/react-table";
import DebouncedInput from "../debounced-input/debounced-input";

type FilterProps = {
  id: string;
  column: Column<any, unknown>;
};

declare module "@tanstack/react-table" {
  //allows us to define custom properties for our columns
  interface ColumnMeta<TData extends RowData, TValue> {
    filterVariant?: "text" | "number";
  }
}

const Filter: React.FunctionComponent<FilterProps> = (props: FilterProps) => {
  const columnFilterValue = props.column.getFilterValue();
  const { filterVariant } = props.column.columnDef.meta ?? {};

  return filterVariant === "number" ? (
    <DebouncedInput
      onChange={(value) => props.column.setFilterValue(value)}
      placeholder={`Search...`}
      type="number"
      value={(columnFilterValue ?? "") as string}
      min={0}
    />
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
