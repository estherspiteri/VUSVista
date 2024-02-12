import React, { useState } from "react";
import styles from "./publication-preview.module.scss";
import { IPublicationPreview } from "../../../models/publication-search/publication-search.model";
import Icon from "../../../atoms/icon/icon";

type PublicationPreviewProps = {
  data?: IPublicationPreview;
};

const PublicationPreview: React.FunctionComponent<PublicationPreviewProps> = (
  props: PublicationPreviewProps
) => {
  const [isAdditionalInfoVisible, setIsAdditionalInfoVisible] = useState(false);

  return (
    <div
      className={`${styles["publication-preview-container"]} ${
        isAdditionalInfoVisible ? styles["visible-additional-info"] : ""
      }`}
    >
      <div
        className={styles.header}
        onClick={() => setIsAdditionalInfoVisible(!isAdditionalInfoVisible)}
      >
        <div className={styles.pmid}>{props.data?.pmid}</div>
        <div className={styles.title}>{props.data?.title}</div>
        <div className={styles["icon-wrapper"]}>
          <div
            className={styles.icon}
            onClick={() =>
              openInNewWindow(
                `https://pubmed.ncbi.nlm.nih.gov/${props.data.pmid}/`
              )
            }
          >
            <Icon name="document" />
          </div>
        </div>
      </div>
      <div className={styles["additional-info"]}>
        <div className={styles["additional-info-content"]}>
          {props.data?.isSupplementaryMaterialMatch && (
            <div className={styles["supplementary-material-match"]}>
              <div className={styles.icon}>
                <Icon name="warning" fill="#fff" width={16} height={16} />
              </div>
              <span>Supplementary material match</span>
            </div>
          )}

          <div className={styles.info}>
            <span className={styles["info-type"]}>Date:</span>
            <span>
              {props.data?.date.getDate()}&nbsp;
              {getMonthString(props.data?.date.getMonth())}
              &nbsp;
              {props.data?.date.getFullYear()}
            </span>
          </div>
          <div className={styles.info}>
            <span className={styles["info-type"]}>Authors:</span>
            <span>{props.data?.authors.join(", ")}</span>
          </div>
          <div className={styles.info}>
            <span className={styles["info-type"]}>Journal:</span>
            <span>{props.data?.journal}</span>
          </div>
          <div className={styles.info}>
            <span className={styles["info-type"]}>Abstract:</span>
            <span>{props.data?.abstract}</span>
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

  function openInNewWindow(url: string) {
    const newWindow = window.open(url, "_blank", "noopener,noreferrer");
    if (newWindow) newWindow.opener = null;
  }
};

export default PublicationPreview;
