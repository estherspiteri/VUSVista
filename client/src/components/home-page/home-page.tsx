import React from "react";
import styles from "./home-page.module.scss";
import Logo from "./../../assets/logo.png";

type HomePageProps = {};

const HomePage: React.FunctionComponent<HomePageProps> = (
  props: HomePageProps
) => {
  return (
    <div className={styles["home-page-container"]}>
      <div className={styles.content}>
        <div className={styles["title-container"]}>
          <img src={Logo} className={styles.logo} />
          <span className={styles.title}>VUSVista</span>
        </div>
      </div>
    </div>
  );
};

export default HomePage;
