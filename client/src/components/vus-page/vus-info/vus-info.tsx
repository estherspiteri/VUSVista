import React from "react";
import styles from "./vus-info.module.scss";
import { IVus } from "../../../models/view-vus.model";

type VusInfoProps = {
  vus: IVus;
};

const VusInfo: React.FunctionComponent<VusInfoProps> = (
  props: VusInfoProps
) => {
  return (
    <div className={styles["vus-info-container"]}>
      <div className={styles["vus-info"]}>
        <div className={styles["top-container"]}>
          {/** General information */}
          <div className={styles.information}>
            <div className={styles.info}>
              <div className={styles["info-title"]}>Vus Id:</div>
              {props.vus.id}
            </div>

            {/* <div className={styles.info}>
              <div className={styles["info-title"]}>Genome Version:</div>
              {props.sample.genomeVersion}
            </div>

            <div className={styles.info}>
              <div className={styles["info-title"]}>File uploads:</div>
              <div>
                {props.sample.files.map((f) => {
                  return (
                    <p>
                      {f.filename}&nbsp;
                      {`(${f.dateOfFileUpload.getDate()}/${
                        f.dateOfFileUpload.getMonth() + 1
                      }/${f.dateOfFileUpload.getFullYear()})`}
                    </p>
                  );
                })}
              </div>
            </div> */}
          </div>
        </div>
      </div>
    </div>
  );
};

export default VusInfo;
