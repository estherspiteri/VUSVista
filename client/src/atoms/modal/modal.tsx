import React from "react";
import styles from "./modal.module.scss";

type ModalProps = {
  title?: string;
  children: React.ReactNode;
};

const Modal: React.FunctionComponent<ModalProps> = (props: ModalProps) => {
  return (
    <>
      <div className={styles["modal-overlay"]} />
      <div className={styles["modal-container"]}>
        {props.title && <p className={styles["modal-title"]}>{props.title}</p>}
        <div>{props.children}</div>
      </div>
    </>
  );
};

export default Modal;
