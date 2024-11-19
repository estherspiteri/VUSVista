import React from "react";
import Text from "../../../../atoms/text/text";
import styles from "./debounced-input.module.scss";

type DebouncedInputProps = {
  value: string | number;
  onChange: (value: string | number) => void;
  debounce?: number;
} & Omit<React.InputHTMLAttributes<HTMLInputElement>, "onChange">;

const DebouncedInput: React.FunctionComponent<DebouncedInputProps> = (
  props: DebouncedInputProps
) => {
  const [inputValue, setInputValue] = React.useState(props.value);

  React.useEffect(() => {
    setInputValue(props.value);
  }, [props.value]);

  React.useEffect(() => {
    const timeout = setTimeout(() => {
      props.onChange(inputValue);
    }, props.debounce);

    return () => clearTimeout(timeout);
  }, [inputValue]);

  return (
    <Text
      {...props}
      value={inputValue}
      onChange={(e) => setInputValue(e.target.value)}
      className={styles["debounced-input"]}
      borderRadius={8}
    />
  );
};

DebouncedInput.defaultProps = {
  debounce: 500,
};

export default DebouncedInput;
