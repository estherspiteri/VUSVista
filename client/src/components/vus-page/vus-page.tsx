import React from "react";
import styles from "./vus-page.module.scss";
import { IVus } from "../../models/view-vus.model";
import VusInfo from "./vus-info/vus-info";
import { IAcmgRule } from "../../models/acmg-rule.model";
import { VusService } from "../../services/vus/vus.service";
import Button from "../../atoms/button/button";
import { Link } from "react-router-dom";

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
        <div>
          <div className={styles.title}>Vus Information</div>
          <div className={styles.description}>
            <p>Below you can find information about the selected variant.</p>
          </div>
        </div>
        <Link to={`/review/${props.vus.id}`}>
          <Button text="Classification Review" icon="edit" />
        </Link>
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
