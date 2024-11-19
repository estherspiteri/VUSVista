import React from "react";
import styles from "./error-page.module.scss";
import { ReactComponent as ErrorSvg } from "./error.svg";

type ErrorPageProps = {};

const ErrorPage: React.FunctionComponent<ErrorPageProps> = (
  props: ErrorPageProps
) => {
  return (
    <div className={styles["error-page-container"]}>
      <div className={styles.content}>
        <ErrorSvg />
        <div className={styles.text}>
          <p className={styles.title}>Oops, something's gone wrong!</p>
          <p className={styles.description}>
            Kindly contact the website administrator and explain what instigated
            this error.
          </p>
        </div>
      </div>
    </div>
  );
};

export default ErrorPage;
