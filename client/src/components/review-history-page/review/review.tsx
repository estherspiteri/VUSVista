import React from "react";
import styles from "./review.module.scss";
import { IClassificationReview } from "../../../models/classification-review.model";
import { getMonthString } from "../../../helpers/date-helper";
import { openInNewWindow } from "../../../helpers/open-links";

type ReviewProps = {
  review?: IClassificationReview;
  dateVariantAdded?: Date;
};

const Review: React.FunctionComponent<ReviewProps> = (props: ReviewProps) => {
  return (
    <div className={styles["review-container"]}>
      <div className={styles.date}>
        {props.review && getDate(props.review.dateAdded)}
        {props.dateVariantAdded && getDate(props.dateVariantAdded)}
      </div>

      {props.review && (
        <div className={styles["review-content"]}>
          {/** Classification */}
          <p
            className={`${styles.classification} ${
              styles[
                props.review.classification.toLowerCase().replace("_", "-")
              ]
            }`}
          >
            {props.review.classification.replace("_", " ")}
          </p>

          {/** ACMG Rules */}
          {props.review.acmgRules.length > 0 && (
            <div className={styles.acmg}>
              <span>
                <b>Relevant ACMG Rules: </b>
              </span>
              {props.review.acmgRules.map((r) => {
                const styleClass = r.substring(0, r.length - 1).toLowerCase();

                return (
                  <div
                    className={`${styles["selected-acmg"]} ${
                      styles[`${styleClass}`]
                    }`}
                  >
                    {r}
                  </div>
                );
              })}
            </div>
          )}

          {/** Publications */}
          {props.review.publications.length > 0 && (
            <div className={styles["publications-container"]}>
              <p>
                <b>Relevant Publications:</b>
              </p>
              <div className={styles.publications}>
                {props.review.publications.map((p) => {
                  return (
                    <div
                      className={`${styles.publication} ${
                        p.link ? styles.link : ""
                      }`}
                      onClick={() => p.link && openInNewWindow(p.link)}
                    >
                      <p>{p.title}</p>
                      <p className={styles.doi}>{p.doi}</p>
                    </div>
                  );
                })}
              </div>
            </div>
          )}

          {/** Reason */}
          {props.review.reason && (
            <div className={styles["reason-container"]}>
              <p>
                <b>Reason: </b>
              </p>
              <p className={styles.reason}>{props.review.reason}</p>
            </div>
          )}

          {/** Submitter */}
          <div className={styles.submitter}>
            <div>
              <b>Submitted by:</b>
            </div>
            <div>
              <p>{props.review.scientificMemberName}</p>
              <a
                href={`mailto:${props.review.scientificMemberEmail}`}
                className={styles.email}
              >
                {props.review.scientificMemberEmail}
              </a>
            </div>
          </div>
        </div>
      )}
      {props.dateVariantAdded && (
        <p className={`${styles.classification} ${styles.vus}`}>VUS</p>
      )}
    </div>
  );

  function getDate(date: Date) {
    return (
      <>
        {date.getDate().toString().length > 1
          ? date.getDate()
          : "0" + date.getDate()}
        &nbsp;
        {getMonthString(date.getMonth())}
        &nbsp;
        {date.getFullYear()}
      </>
    );
  }
};

export default Review;
