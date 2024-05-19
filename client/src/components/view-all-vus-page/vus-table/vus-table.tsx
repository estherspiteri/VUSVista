import React from "react";
import styles from "./vus-table.module.scss";
import ViewVus from "../view-vus/view-vus";
import { IVUSSummary } from "../../../models/vus-summary.model";
import {
  ColumnDef,
  ColumnFiltersState,
  Row,
  SortingState,
  flexRender,
  getCoreRowModel,
  getFilteredRowModel,
  getSortedRowModel,
  useReactTable,
} from "@tanstack/react-table";
import Filter from "../../shared/table/filter/filter";
import Icon from "../../../atoms/icons/icon";

type VusTableProps = {
  vusList: IVUSSummary[];
};

const VusTable: React.FunctionComponent<VusTableProps> = (
  props: VusTableProps
) => {
  const [data, setData] = React.useState<IVUSSummary[]>(props.vusList);
  const [columnFilters, setColumnFilters] = React.useState<ColumnFiltersState>(
    []
  );
  const [sorting, setSorting] = React.useState<SortingState>([]);

  //define columns
  const columns = React.useMemo<ColumnDef<IVUSSummary, any>[]>(
    () => [
      {
        accessorKey: "id",
        cell: (info) => info.getValue(),
        header: () => <span>Variant Id</span>,
        meta: {
          filterVariant: "number",
        },
        filterFn: "includesString",
      },
      {
        accessorKey: "chromosome",
        cell: (info) => info.getValue(),
        header: () => <span>Chromosome</span>,
        sortingFn: chromosomeSortFn,
      },
      {
        accessorKey: "chromosomePosition",
        cell: (info) => info.getValue(),
        header: () => <span>Position</span>,
        meta: {
          filterVariant: "number",
        },
        filterFn: "includesString",
        sortingFn: "alphanumeric",
      },
      {
        accessorKey: "gene",
        cell: (info) => info.getValue(),
        header: () => <span>Gene</span>,
      },
      {
        accessorKey: "refAllele",
        cell: (info) => info.getValue(),
        header: () => <span>Reference</span>,
      },
      {
        accessorKey: "altAllele",
        cell: (info) => info.getValue(),
        header: () => <span>Alternate</span>,
      },
      {
        accessorKey: "rsid",
        cell: (info) => info.getValue(),
        header: () => <span>RSID</span>,
        sortingFn: rsidSortFn,
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
    getFilteredRowModel: getFilteredRowModel(), //client side filtering
    onColumnFiltersChange: setColumnFilters,
    getSortedRowModel: getSortedRowModel(), //client-side sorting
    onSortingChange: setSorting, //optionally control sorting state in your own scope for easy access
  });

  return (
    <div className={styles["vus-table-container"]}>
      <table className={styles.table}>
        <thead>
          {table.getHeaderGroups().map((headerGroup) => (
            <tr key={headerGroup.id} className={styles.header}>
              {headerGroup.headers.map((header) => {
                return (
                  <th key={header.id} className={styles["header-content"]}>
                    {header.isPlaceholder ? null : (
                      <>
                        <div className={styles.title}>
                          {flexRender(
                            header.column.columnDef.header,
                            header.getContext()
                          )}
                          <div
                            className={styles["sort-icons"]}
                            onClick={header.column.getToggleSortingHandler()}
                            title={
                              header.column.getCanSort()
                                ? header.column.getNextSortingOrder() === "asc"
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
                              desc: <Icon name="desc" width={16} height={16} />,
                              false: (
                                <Icon name="sort" width={16} height={16} />
                              ),
                            }[header.column.getIsSorted() as string] ?? null}
                          </div>
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
            <ViewVus vusRow={row} isColoured={index % 2 === 0} />
          ))}
        </tbody>
      </table>
    </div>
  );

  function chromosomeSortFn(
    rowA: Row<IVUSSummary>,
    rowB: Row<IVUSSummary>,
    columnId: string
  ): number {
    const a: string = rowA.renderValue(columnId);
    const b: string = rowB.renderValue(columnId);

    var numCheckRegex = /^\d+$/;

    const aIsNumber = a.match(numCheckRegex);
    const bIsNumber = b.match(numCheckRegex);

    if (aIsNumber && bIsNumber) {
      return a > b ? 1 : -1;
    } else if (aIsNumber) {
      return -1;
    } else if (bIsNumber) {
      return 1;
    } else {
      return a > b ? 1 : -1;
    }
  }

  function rsidSortFn(
    rowA: Row<IVUSSummary>,
    rowB: Row<IVUSSummary>,
    columnId: string
  ): number {
    //ignore the 'rs' prefix
    const a: string = rowA.renderValue(columnId) ?? "".split("rs")[1];
    const b: string = rowB.renderValue(columnId) ?? "".split("rs")[1];

    return a > b ? 1 : -1;
  }
};

export default VusTable;
