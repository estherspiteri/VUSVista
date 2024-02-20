import React from "react";
import styles from "./publication-view-page.module.scss";
import PublicationPreview from "./publication-preview/publication-preview";
import { IPublicationPreview } from "../../models/publication-search.model";

type PublicationViewPageProps = {
  variantId: string;
  publications?: IPublicationPreview[];
};

const PublicationViewPage: React.FunctionComponent<PublicationViewPageProps> = (
  props: PublicationViewPageProps
) => {
  return (
    <div className={styles["publication-view-container"]}>
      <div className={styles["title-container"]}>
        <div className={styles.title}>Publications</div>
        <div className={styles.description}>
          <p>Look up publications for a particular VUS using its Variant Id.</p>
        </div>
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
                &nbsp;found for Variant with &nbsp;
                <span className={styles.colour}>Id {props.variantId}</span>
              </span>
              <div className={styles["publication-previews"]}>
                <div className={styles.header}>
                  <span className={styles.pmid}>PMID</span>
                  <span>Title</span>
                </div>

                <div className={styles["publication-preview-contents"]}>
                  {props.publications.map((publication) => (
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
};

export default PublicationViewPage;
