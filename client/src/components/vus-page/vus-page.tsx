import React from "react";
import styles from "./vus-page.module.scss";
import { IVus } from "../../models/view-vus.model";
import VusInfo from "./vus-info/vus-info";
import { IAcmgRule } from "../../models/acmg-rule.model";
import { VusService } from "../../services/vus/vus.service";
import { Link } from "react-router-dom";
import Icon from "../../atoms/icons/icon";

type VusPageProps = {
  vus: IVus;
  acmgRules: IAcmgRule[];
  vusService?: VusService;
};

const VusPage: React.FunctionComponent<VusPageProps> = (
  props: VusPageProps
) => {
  return (
    <div className={styles["vus-page-container"]}>
      <div className={styles["title-container"]}>
        <div className={styles["title-section"]}>
          <div className={styles.title}>Variant Information</div>
          <div className={styles["option-container"]}>
            <div className={styles["option-btn"]}>
              <span>Variant Options</span>
              <Icon name="options" width={12} height={12} />
            </div>
            <div className={styles.options}>
              <Link to={`/review/${props.vus.id}`} className={styles.option}>
                New Classification Review
              </Link>
              <Link
                to={`/review-history/${props.vus.id}`}
                className={styles.option}
              >
                Classification Review History
              </Link>
              {props.vus.numOfPublications > 0 && (
                <Link
                  to={`/publication-view/${props.vus.id}`}
                  className={styles.option}
                >
                  View Publications
                </Link>
              )}
            </div>
          </div>
        </div>
        <div className={styles.description}>
          <p>Below you can find information about the selected variant.</p>
        </div>
      </div>

      <VusInfo
        vus={props.vus}
        acmgRules={props.acmgRules}
        vusService={props.vusService}
      />
    </div>
  );
};

export default VusPage;
