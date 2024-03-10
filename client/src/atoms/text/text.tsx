import React, { useState } from "react";
import styles from "./text.module.scss";

type TextProps = {
  errorMsg?: string;
  validationCallback?: (val: string) => boolean;
} & React.InputHTMLAttributes<HTMLInputElement>;

const Text: React.FunctionComponent<TextProps> = (props: TextProps) => {
  const [isValid, setIsValid] = useState(props.errorMsg.length > 0);

  const { validationCallback, ...remainingProps } = props;

  return (
    <div className={`${styles.text} ${props.className}`}>
      <div className={styles["input-container"]}>
        <input
          onBlur={(e) => {
            if (props.validationCallback) {
              setIsValid(props.validationCallback(e.currentTarget.value));
            }
          }}
          {...remainingProps}
        />
      </div>
      {props.errorMsg?.length > 0 && (
        <p className={styles.error}>{props.errorMsg}</p>
      )}
    </div>
  );
};

export default Text;
