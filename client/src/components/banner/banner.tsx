import React, { FunctionComponent, PropsWithChildren, useState } from "react";
import styles from "./banner.module.scss";

type BannerProps = { isClosable: boolean } & PropsWithChildren;

export const Banner: FunctionComponent<BannerProps> = (props: BannerProps) => {
  const [isOpen, setIsOpen] = useState(true);

  return (
    <div
      className={`${styles["banner-container"]} ${isOpen ? "" : styles.closed}`}
    >
      <div className={styles.content}>{props.children}</div>
      <div onClick={() => setIsOpen(false)} className={styles["close-icon"]}>
        <svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
          <path
            d="m16 16-4-4m0 0L8 8m4 4 4-4m-4 4-4 4"
            stroke="#fff"
            stroke-width="2"
            stroke-linecap="round"
            stroke-linejoin="round"
          />
        </svg>
      </div>
    </div>
  );
};

export default Banner;
