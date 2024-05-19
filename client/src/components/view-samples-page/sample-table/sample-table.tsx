import React from "react";
import styles from "./sample-table.module.scss";
import ViewSample from "../view-sample/view-sample";
import { ISampleSummary } from "../../../models/view-samples.model";
import {
  useReactTable,
  getCoreRowModel,
  flexRender,
  getFilteredRowModel,
  ColumnDef,
  ColumnFiltersState,
  SortingState,
  getSortedRowModel,
} from "@tanstack/react-table";
import Filter from "../../shared/table/filter/filter";
import Icon from "../../../atoms/icons/icon";

type SampleTableProps = {
  sampleList: ISampleSummary[];
};

const SampleTable: React.FunctionComponent<SampleTableProps> = (
  props: SampleTableProps
) => {
  const [data, setData] = React.useState<ISampleSummary[]>(props.sampleList);
  const [columnFilters, setColumnFilters] = React.useState<ColumnFiltersState>(
    []
  );
  const [sorting, setSorting] = React.useState<SortingState>([]);

  //define columns
  const columns = React.useMemo<ColumnDef<ISampleSummary, any>[]>(
    () => [
      {
        accessorKey: "sampleId",
        cell: (info) => info.getValue(),
        header: () => <span>Sample Id</span>,
        enableSorting: false,
      },
      {
        accessorKey: "numOfVariants",
        cell: (info) => info.getValue(),
        header: () => <span>Number of Variants</span>,
        meta: {
          filterVariant: "number",
        },
        filterFn: "weakEquals",
      },
    ],
    []
  );

  const table = useReactTable({
    columns,
    data,
    filterFns: {},
    state: {
      columnFilters,
      sorting,
    },
    getCoreRowModel: getCoreRowModel(),
    onColumnFiltersChange: setColumnFilters,
    getFilteredRowModel: getFilteredRowModel(), //client side filtering
    getSortedRowModel: getSortedRowModel(), //client-side sorting
    onSortingChange: setSorting, //optionally control sorting state in your own scope for easy access
  });

  return (
    <div className={styles["sample-table-container"]}>
      <table>
        <thead>
          {table.getHeaderGroups().map((headerGroup) => (
            <tr key={headerGroup.id} className={styles.header}>
              {headerGroup.headers.map((header) => {
                return (
                  <th
                    key={header.id}
                    className={`${styles["header-content"]} ${
                      header.id === "sampleId" ? styles.id : styles.variants
                    }`}
                  >
                    {header.isPlaceholder ? null : (
                      <>
                        <div className={styles.title}>
                          {flexRender(
                            header.column.columnDef.header,
                            header.getContext()
                          )}
                          {header.column.getCanSort() && (
                            <div
                              className={styles["sort-icons"]}
                              onClick={header.column.getToggleSortingHandler()}
                              title={
                                header.column.getCanSort()
                                  ? header.column.getNextSortingOrder() ===
                                    "asc"
                                    ? "Sort ascending"
                                    : header.column.getNextSortingOrder() ===
                                      "desc"
                                    ? "Sort descending"
                                    : "Clear sort"
                                  : undefined
                              }
                            >
                              {{
                                asc: <Icon name="asc" width={16} height={16} />,
                                desc: (
                                  <Icon name="desc" width={16} height={16} />
                                ),
                                false: (
                                  <Icon name="sort" width={16} height={16} />
                                ),
                              }[header.column.getIsSorted() as string] ?? null}
                            </div>
                          )}
                        </div>
                        {header.column.getCanFilter() ? (
                          <div>
                            <Filter column={header.column} id={header.id} />
                          </div>
                        ) : null}
                      </>
                    )}
                  </th>
                );
              })}
            </tr>
          ))}
        </thead>
        <tbody>
          {table.getRowModel().rows.map((row, index) => (
            <ViewSample sampleRow={row} isColoured={index % 2 === 0} />
          ))}
        </tbody>
      </table>
    </div>
  );
};

export default SampleTable;
