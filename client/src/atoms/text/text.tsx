import React, { CSSProperties, useEffect, useState } from "react";
import styles from "./text.module.scss";

type TextProps = {
  value?: string | number;
  debounce?: number;
  errorMsg?: string;
  borderRadius?: number;
  onValueChange?: (value: string | number) => void;
  validationCallback?: (val: string) => boolean;
} & React.InputHTMLAttributes<HTMLInputElement>;

const Text: React.FunctionComponent<TextProps> = (props: TextProps) => {
  const [isValid, setIsValid] = useState(props.errorMsg?.length > 0);
  const [value, setValue] = useState(props.value);

  const { validationCallback, onChange, ...remainingProps } = props;

  //throttle
  useEffect(() => {
    const timeout = setTimeout(() => {
      if (props.onValueChange && props.debounce) props.onValueChange(value);
    }, props.debounce);

    return () => clearTimeout(timeout);
  }, [value]);

  return (
    <div
      className={`${styles.text} ${props.className}`}
      style={
        { "--border-radius": `${props.borderRadius ?? 4}px` } as CSSProperties
      }
    >
      <div className={styles["input-container"]}>
        <input
          value={value}
          onBlur={(e) => {
            if (props.validationCallback) {
              setIsValid(props.validationCallback(e.currentTarget.value));
            }
          }}
          step={1}
          onChange={(e) => {
            setValue(e.target.value);

            if (props.type === "number" && props.validationCallback) {
              setIsValid(props.validationCallback(e.currentTarget.value));
            }

            onChange(e);
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
