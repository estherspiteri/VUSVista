import React, { FunctionComponent, PropsWithChildren, useState } from "react";
import styles from "./banner.module.scss";
import Icon from "../icons/icon";

type BannerProps = {
  isClosable: boolean;
  isGreen?: boolean;
  onCloseCallback?: () => void;
} & PropsWithChildren;

export const Banner: FunctionComponent<BannerProps> = (props: BannerProps) => {
  const [isOpen, setIsOpen] = useState(true);

  return (
    <div
      className={`${styles["banner-container"]} ${
        isOpen ? "" : styles.closed
      } ${props.isGreen ? styles.green : ""}`}
    >
      <div className={styles.content}>{props.children}</div>
      {/* <div onClick={() => setIsOpen(false)} className={styles["close-icon"]}> */}
      <Icon
        name="close"
        stroke="#fff"
        className={styles["close-icon"]}
        onClick={() => {
          setIsOpen(false);
          props.onCloseCallback && props.onCloseCallback();
        }}
      />
      {/* </div> */}
    </div>
  );
};

export default Banner;
