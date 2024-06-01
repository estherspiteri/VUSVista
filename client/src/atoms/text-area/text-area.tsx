import React, { CSSProperties, useState } from "react";
import styles from "./text-area.module.scss";

type TextAreaProps = {
  value?: string | number;
  borderRadius?: number;
} & React.TextareaHTMLAttributes<HTMLTextAreaElement>;

const TextArea: React.FunctionComponent<TextAreaProps> = (
  props: TextAreaProps
) => {
  const [value, setValue] = useState(props.value);

  const { borderRadius, onChange, ...remainingProps } = props;

  return (
    <div
      className={`${styles["text-area-container"]} ${props.className}`}
      style={{ "--border-radius": `${borderRadius ?? 4}px` } as CSSProperties}
    >
      <textarea
        value={value}
        onChange={(e) => {
          setValue(e.target.value);
          onChange(e);
        }}
        {...remainingProps}
      />
    </div>
  );
};

export default TextArea;
