import React from "react";
import styles from "./review-history-page.module.scss";
import { IVUSSummary } from "../../models/vus-summary.model";
import VariantSummary from "../shared/variant-summary/variant-summary";
import { IClassificationReview } from "../../models/classification-review.model";
import Review from "./review/review";

type ReviewHistoryPageProps = {
  variantSummary: IVUSSummary;
  reviews: IClassificationReview[];
};

const ReviewHistoryPage: React.FunctionComponent<ReviewHistoryPageProps> = (
  props: ReviewHistoryPageProps
) => {
  return (
    <div className={styles["review-history-container"]}>
      <div className={styles["title-container"]}>
        <div className={styles.title}>Classification Review History</div>
        <div className={styles.description}>
          <p>
            This page contains the below variant's creation date and any submitted
            Classification Reviews that might have altered this
            variant's classification.
          </p>
        </div>
      </div>
      <div className={styles["variant-info"]}>
        <div className={styles["variant-summary"]}>
          <VariantSummary variant={props.variantSummary} />
        </div>

        <div className={styles.reviews}>
          {props.reviews.map((r) => {
            return <Review review={r} />;
          })}
        </div>
      </div>
    </div>
  );
};

export default ReviewHistoryPage;
