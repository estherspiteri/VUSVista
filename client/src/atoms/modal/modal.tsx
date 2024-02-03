import React, {
  ForwardRefRenderFunction,
  PropsWithChildren,
  forwardRef,
  useImperativeHandle,
  useState,
} from "react";
import styles from "./modal.module.scss";

export interface ModalRef {
  closeModal: () => void;
}

type ModalProps = {
  title?: string;
};

const ModalRef: ForwardRefRenderFunction<
  ModalRef,
  PropsWithChildren<ModalProps>
> = (props, ref) => {
  const [isModalOpen, setIsModalOpen] = useState(true);

  //expose internal methods
  useImperativeHandle(ref, () => ({ closeModal: closeModal }));

  return (
    <div style={{ display: isModalOpen ? "block" : "none" }}>
      <div className={styles["modal-overlay"]} />
      <div className={styles["modal-container"]}>
        {props.title && <p className={styles["modal-title"]}>{props.title}</p>}
        <div>{props.children}</div>
      </div>
    </div>
  );

  function closeModal() {
    setIsModalOpen(false);
  }
};

export const Modal = forwardRef(ModalRef);

export default Modal;
