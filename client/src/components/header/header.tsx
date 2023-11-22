import React from "react";
import styles from "./header.module.scss";
import { Link } from "react-router-dom";

type HeaderProps = {};

const Header: React.FunctionComponent<HeaderProps> = (props: HeaderProps) => {
  return (
    <div className={styles["header-container"]}>
      <div className={styles["header-content"]}>
        <Link to={""} className={styles.option}>
          VUS file upload
        </Link>
        <Link to={"/view-vus"} className={styles.option}>
          View all VUS
        </Link>
        <Link to={"/publication-search"} className={styles.option}>
          Publication search
        </Link>
      </div>
    </div>
  );
};

export default Header;
