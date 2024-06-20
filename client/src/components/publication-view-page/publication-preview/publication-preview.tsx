import React, { useEffect, useRef, useState } from "react";
import styles from "./publication-preview.module.scss";
import Icon from "../../../atoms/icons/icon";
import { IPublicationPreview } from "../../../models/publication-view.model";
import { openInNewWindow } from "../../../helpers/open-links";
import { getMonthString } from "../../../helpers/date-helper";

type PublicationPreviewProps = {
  data?: IPublicationPreview;
};

const PublicationPreview: React.FunctionComponent<PublicationPreviewProps> = (
  props: PublicationPreviewProps
) => {
  const ref = useRef<HTMLDivElement>(null);
  const [isAdditionalInfoVisible, setIsAdditionalInfoVisible] = useState(false);

  //close additional info on outside click
  useEffect(() => {
    function handleClickOutside(event) {
      if (ref.current && !ref.current.contains(event.target)) {
        setIsAdditionalInfoVisible(false);
      }
    }
    // Bind the event listener
    document.addEventListener("mousedown", handleClickOutside);

    return () => {
      // Unbind the event listener on clean up
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [ref]);

  return (
    <div
      ref={ref}
      className={`${styles["publication-preview-container"]} ${
        isAdditionalInfoVisible ? styles["visible-additional-info"] : ""
      }`}
    >
      <div
        className={`${styles.header} ${
          isAdditionalInfoVisible ? styles["header-with-add-info"] : ""
        }`}
        onClick={() => setIsAdditionalInfoVisible(!isAdditionalInfoVisible)}
      >
        <div className={styles.title}>{props.data?.title}</div>
        {props.data.isAddedManually && (
          <div className={styles["profile-icon-wrapper"]}>
            <Icon name="profile" fill="#008080" />
          </div>
        )}
        <div className={styles["external-link-icon-wrapper"]}>
          <div
            className={styles.icon}
            onClick={() => openInNewWindow(props.data.link)}
          >
            <Icon
              name="external-link"
              className={styles["external-link-icon"]}
            />
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
          {/** Publication Information */}
          <div className={styles.info}>
            <span className={styles["info-type"]}>DOI:</span>
            <span>{props.data.doi}</span>
          </div>
          {props.data.pmid && (
            <div className={styles.info}>
              <span className={styles["info-type"]}>PMID:</span>
              <span>{props.data.pmid}</span>
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
          {props.data.authors && (
            <div className={styles.info}>
              <span className={styles["info-type"]}>Authors:</span>
              <span>
                {props.data?.authors.slice(0, 10)?.join(", ")}
                {props.data.authors.length > 9 ? ", ..." : ""}
              </span>
            </div>
          )}
          {props.data.journal && (
            <div className={styles.info}>
              <span className={styles["info-type"]}>Journal:</span>
              <span>{props.data?.journal}</span>
            </div>
          )}
          {props.data.abstract && (
            <div className={styles.info}>
              <span className={styles["info-type"]}>Abstract:</span>
              <p dangerouslySetInnerHTML={{ __html: props.data?.abstract }} />
            </div>
          )}
          <div className={styles["publication-url"]}>
            {props.data.link && (
              <div className={styles.info}>
                <span className={styles["info-type"]}>URL:</span>
                <a href={props.data?.link}>{props.data?.link}</a>
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default PublicationPreview;
