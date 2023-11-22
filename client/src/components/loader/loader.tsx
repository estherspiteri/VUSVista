import React from "react";
import styles from "./loader.module.scss";

type LoaderProps = {};

const Loader: React.FunctionComponent<LoaderProps> = (props: LoaderProps) => {
  return <div className={styles["lds-dual-ring"]}></div>;
};

export default Loader;
