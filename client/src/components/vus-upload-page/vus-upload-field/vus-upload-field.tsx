import React, { PropsWithChildren, useEffect, useRef, useState } from "react";
import styles from "./vus-upload-field.module.scss";
import Icon from "../../../atoms/icons/icon";

type VusUploadFieldProps = PropsWithChildren<{
  title: string;
  showCheckMark?: boolean;
  isOpenByDefault?: boolean;
}>;

const VusUploadField: React.FunctionComponent<VusUploadFieldProps> = (
  props: VusUploadFieldProps
) => {
  const ref = useRef<HTMLDivElement>(null);

  const [isContentVisible, setIsContentVisible] = useState(
    props.isOpenByDefault ?? false
  );

  //close content on outside click
  useEffect(() => {
    function handleClickOutside(event) {
      if (ref.current && !ref.current.contains(event.target)) {
        setIsContentVisible(false);
      }
    }
    // Bind the event listener
    document.addEventListener("mousedown", handleClickOutside);

    return () => {
      // Unbind the event listener on clean up
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [ref]);

  return (
    <div ref={ref} className={styles["vus-upload-field-container"]}>
      <div
        className={`${styles.header} ${isContentVisible ? styles.open : ""}`}
        onClick={() => {
          setIsContentVisible(!isContentVisible);
        }}
      >
        <div className={styles["header-left"]}>
          <Icon name="chev-right" className={styles["chev-right"]} />
          <span>{props.title}</span>
        </div>
        {props.showCheckMark && <Icon name="checkmark" />}
      </div>
      <div
        className={`${styles.content} ${
          isContentVisible ? styles["content-visible"] : ""
        }`}
      >
        <div className={styles.children}>{props.children}</div>
      </div>
    </div>
  );
};

export default VusUploadField;
