import React, { useState } from "react";
import WarningIcon from "./../../../src/warning.svg";
import styles from "./literature-preview.module.scss";

type LiteraturePreviewProps = {
  pmid: number;
  title: string;
  date: Date;
  authors: string[];
  journal: string;
  abstract: string;
  isSupplementaryMaterialMatch: boolean;
};

const LiteraturePreview: React.FunctionComponent<LiteraturePreviewProps> = (
  props: LiteraturePreviewProps
) => {
  const [isAdditionalInfoVisible, setIsAdditionalInfoVisible] = useState(false);

  return (
    <div className={styles["literature-preview-component"]}>
      <div
        className={styles.header}
        onClick={() => setIsAdditionalInfoVisible(!isAdditionalInfoVisible)}
      >
        <div className={styles.pmid}>{props.pmid}</div>
        <div className={styles.title}>{props.title}</div>
      </div>
      <div
        className={`${styles["additional-info"]} ${
          isAdditionalInfoVisible ? styles["visible-additional-info"] : ""
        }`}
      >
        <div className={styles["additional-info-content"]}>
          {props.isSupplementaryMaterialMatch && (
            <div className={styles["supplementary-material-match"]}>
              <div className={styles.icon}>
                <svg fill="white" viewBox="0 0 478.125 478.125">
                  <circle cx="239.904" cy="314.721" r="35.878" />
                  <path d="M256.657,127.525h-31.9c-10.557,0-19.125,8.645-19.125,19.125v101.975c0,10.48,8.645,19.125,19.125,19.125h31.9 c10.48,0,19.125-8.645,19.125-19.125V146.65C275.782,136.17,267.138,127.525,256.657,127.525z" />
                  <path d="M239.062,0C106.947,0,0,106.947,0,239.062s106.947,239.062,239.062,239.062c132.115,0,239.062-106.947,239.062-239.062 S371.178,0,239.062,0z M239.292,409.734c-94.171,0-170.595-76.348-170.595-170.596c0-94.248,76.347-170.595,170.595-170.595 s170.595,76.347,170.595,170.595C409.887,333.387,333.464,409.734,239.292,409.734z" />
                </svg>
              </div>
              <span>Supplementary material match</span>
            </div>
          )}

          <div className={styles.info}>
            <span className={styles["info-type"]}>Date:</span>
            <span>
              {props.date.getDate()} {getMonthString(props.date.getMonth())}
              &nbsp;
              {props.date.getFullYear()}
            </span>
          </div>
          <div className={styles.info}>
            <span className={styles["info-type"]}>Authors:</span>
            <span>{props.authors.join(", ")}</span>
          </div>
          <div className={styles.info}>
            <span className={styles["info-type"]}>Journal:</span>
            <span>{props.journal}</span>
          </div>
          <div className={styles.info}>
            <span className={styles["info-type"]}>Abstract:</span>
            <span>{props.abstract}</span>
          </div>
        </div>
      </div>
    </div>
  );

  function getMonthString(month: number) {
    const months = [
      "Jan",
      "Feb",
      "Mar",
      "Apr",
      "May",
      "Jun",
      "Jul",
      "Aug",
      "Sep",
      "Oct",
      "Nov",
      "Dec",
    ];

    return months[month];
  }
};

export default LiteraturePreview;
