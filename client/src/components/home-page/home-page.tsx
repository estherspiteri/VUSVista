import React from "react";
import styles from "./home-page.module.scss";
import Logo from "./../../assets/logo.png";
import { IHomepageData } from "../../models/homepage.model";
import VusTable from "../view-all-vus-page/vus-table/vus-table";

type HomePageProps = { data: IHomepageData };

const HomePage: React.FunctionComponent<HomePageProps> = (
  props: HomePageProps
) => {
  return (
    <>
      <div className={styles.banner}>
        <div className={styles["title-container"]}>
          <img src={Logo} className={styles.logo} />
          <span className={styles.title}>VUSVista</span>
        </div>
      </div>
      <div className={styles["home-page-container"]}>
        <div className={styles.content}>
          <div className={styles.updates}>
            <p className={styles["last-update"]}>
              <span className={styles["last-update-desc"]}>
                Clinvar last auto update:
              </span>
              &nbsp;
              {props.data.lastClinvarUpdateDate}
            </p>
            <p className={styles["last-update"]}>
              <span className={styles["last-update-desc"]}>
                Publications last auto update:
              </span>
              &nbsp;
              {props.data.lastPubUpdateDate}
            </p>
          </div>

          {props.data.vusList && (
            <div>
              <p className={styles["info-title"]}>Latest uploaded variants</p>
              <VusTable
                vusList={props.data.vusList}
                showFilters={false}
                showSort={false}
                showExtraInfoColumns={false}
              />
            </div>
          )}
        </div>
      </div>
    </>
  );
};

export default HomePage;
