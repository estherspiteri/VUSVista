import React, { useState } from "react";
import styles from "./header.module.scss";
import { Link } from "react-router-dom";

type HeaderProps = {};

const Header: React.FunctionComponent<HeaderProps> = (props: HeaderProps) => {
  return (
    <div className={styles["header-container"]}>
      <div className={styles["header-content"]}>
        {/** VUS */}
        <div
          className={`${styles["option-container"]} ${styles["vus-option-container"]}`}
        >
          <div className={styles["option-btn"]}>VUS</div>
          <div className={styles.options}>
            <Link to={"upload-vus"} className={styles.option}>
              File Upload
            </Link>
            <Link to={"/view-vus"} className={styles.option}>
              View All
            </Link>
          </div>
        </div>

        {/** Publications */}
        <div
          className={`${styles["option-container"]} ${styles["publications-option-container"]}`}
        >
          <div className={styles["option-btn"]}>Publications</div>
          <div className={styles.options}>
            <Link to={"/publication-search"} className={styles.option}>
              Search
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
};

export default Header;
