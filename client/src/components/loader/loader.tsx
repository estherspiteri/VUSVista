import React, { CSSProperties } from "react";
import styles from "./loader.module.scss";

type LoaderProps = { width?: number; thickness?: number };

const Loader: React.FunctionComponent<LoaderProps> = (props: LoaderProps) => {
  return (
    <div
      style={
        {
          "--width": `${props.width ?? 64}px`,
          "--thick": `${props.thickness ?? 6}px`,
        } as CSSProperties
      }
      className={styles["lds-dual-ring"]}
    ></div>
  );
};

export default Loader;
