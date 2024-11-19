import React, {
  ForwardRefRenderFunction,
  PropsWithChildren,
  forwardRef,
  useImperativeHandle,
  useState,
} from "react";
import styles from "./modal.module.scss";
import Icon from "../icons/icon";

export interface ModalRef {
  closeModal: () => void;
}

type ModalProps = {
  isClosable?: boolean;
  title?: string;
  modalContainerStyle?: string;
  onCloseIconClickCallback?: () => void;
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
      <div
        className={`${styles["modal-container"]} ${props.modalContainerStyle}`}
      >
        {props.title && (
          <div
            className={`${styles["modal-title"]} ${
              props.isClosable ? styles.closable : ""
            }`}
          >
            <span>{props.title}</span>
            {props.isClosable && (
              <Icon
                name="close"
                onClick={props.onCloseIconClickCallback}
                width={32}
                height={32}
                className={styles["close-icon"]}
              />
            )}
          </div>
        )}
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
