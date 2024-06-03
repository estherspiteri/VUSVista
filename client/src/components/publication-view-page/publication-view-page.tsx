import React, { useState } from "react";
import styles from "./publication-view-page.module.scss";
import PublicationPreview from "./publication-preview/publication-preview";
import { IPublicationPreview } from "../../models/publication-view.model";
import VariantSummary from "../shared/variant-summary/variant-summary";
import { IVUSSummary } from "../../models/vus-summary.model";
import Text from "../../atoms/text/text";
import DebouncedInput from "../shared/table/debounced-input/debounced-input";

type PublicationViewPageProps = {
  description?: string;
  variantId: string;
  variant: IVUSSummary;
  publications?: IPublicationPreview[];
};

const PublicationViewPage: React.FunctionComponent<PublicationViewPageProps> = (
  props: PublicationViewPageProps
) => {
  const [filterValue, setFilterValue] = useState("");

  return (
    <div className={styles["publication-view-container"]}>
      <div className={styles["title-container"]}>
        <div className={styles.title}>Publications</div>
        {props.description && (
          <div
            className={styles.description}
            dangerouslySetInnerHTML={{ __html: props.description }}
          />
        )}
      </div>

      {props.publications && (
        <>
          {props.publications && (
            <div className={styles["publications-previews-container"]}>
              <span className={styles.status}>
                <span className={styles.colour}>
                  {props.publications.length}
                </span>
                &nbsp;
                {props.publications.length === 1
                  ? "publication"
                  : "publications"}
                &nbsp;found for the below variant
              </span>
              <div className={styles["variant-summary"]}>
                <VariantSummary variant={props.variant} />
              </div>
              <div className={styles["publication-previews"]}>
                <div className={styles.header}>
                  <span className={styles["header-title"]}>
                    Publication Titles
                  </span>
                  <DebouncedInput
                    onChange={(val) => setFilterValue(val.toString())}
                    placeholder={`Search publication titles...`}
                    type="text"
                    value={filterValue}
                    className={styles.input}
                  />
                </div>

                <div className={styles["publication-preview-contents"]}>
                  {(filterValue.length > 0
                    ? props.publications.filter(
                        (p) =>
                          containsFilterValue(p.title) ||
                          containsFilterValue(p.doi) ||
                          containsFilterValue(p.abstract) ||
                          (p.authors &&
                            containsFilterValue(p.authors.join(" "))) ||
                          containsFilterValue(p.journal) ||
                          containsFilterValue(p.pmid)
                      )
                    : props.publications
                  ).map((publication) => (
                    <PublicationPreview data={publication} />
                  ))}
                </div>
              </div>
            </div>
          )}
        </>
      )}
    </div>
  );

  function containsFilterValue(val: string): boolean {
    if (val && val.length > 0) {
      return val?.toLowerCase().includes(filterValue?.toLowerCase());
    }
    return false;
  }
};

export default PublicationViewPage;
