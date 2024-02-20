import React, { useState } from "react";
import styles from "./header.module.scss";
import { Link } from "react-router-dom";

type HeaderProps = {};

const Header: React.FunctionComponent<HeaderProps> = (props: HeaderProps) => {
  return (
    <div className={styles["header-container"]}>
      <div className={styles["header-content"]}>
        {/** VUS */}
        <div className={styles["option-container"]}>
          <div className={styles["option-btn"]}>VUS</div>
          <div className={styles.options}>
            <Link to={"/view-vus"} className={styles.option}>
              View All
            </Link>
          </div>
        </div>

        {/** Samples */}
        <div className={styles["option-container"]}>
          <div className={styles["option-btn"]}>Samples</div>
          <div className={styles.options}>
            <Link to={"upload-vus"} className={styles.option}>
              File Upload
            </Link>
            <Link to={"/view-samples"} className={styles.option}>
              View All
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
};

export default Header;
