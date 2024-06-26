import React from "react";
import styles from "./home-page.module.scss";

type HomePageProps = {};

const HomePage: React.FunctionComponent<HomePageProps> = (
  props: HomePageProps
) => {
  return (
    <div className={styles["home-page-container"]}>
      <div className={styles.content}>
        This is the home page
      </div>
    </div>
  );
};

export default HomePage;
