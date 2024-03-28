import React from "react";
import styles from "./button.module.scss";
import Icon, { IconName } from "../icons/icon";

export type ButtonType = "default";

type ButtonProps = {
  text: string;
  buttonType?: ButtonType;
  icon?: IconName;
} & React.ButtonHTMLAttributes<HTMLButtonElement>;

const Button: React.FunctionComponent<ButtonProps> = (props: ButtonProps) => {
  const { className, ...remainingProps } = props;

  return (
    <div className={props.className}>
      <button className={styles.button} {...remainingProps}>
        <div className={styles["btn-content"]}>
          {props.text}
          {props.icon && <Icon name={props.icon} width={16} height={16} />}
        </div>
      </button>
    </div>
  );
};

Button.defaultProps = { buttonType: "default" };

export default Button;
