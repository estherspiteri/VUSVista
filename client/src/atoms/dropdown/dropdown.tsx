import React, { CSSProperties, useEffect, useRef, useState } from "react";
import styles from "./dropdown.module.scss";
import Icon from "../icons/icon";

type DropdownElt = { elt: {}; displayElt: JSX.Element };

type DropdownProps = {
  isDisabled?: boolean;
  inputPlaceholder?: string;
  list?: { elt: {}; displayElt: JSX.Element }[];
  openOnClick?: boolean;
  borderRadius?: number;
  dropdownContentContainerClassname?: string;
  onEltClickCallback?: (elt: {}) => void;
  onInputFocusCallback?: (val: string) => void;
  onInputChangeCallback?: (val: string) => void;
};

const Dropdown: React.FunctionComponent<DropdownProps> = (
  props: DropdownProps
) => {
  const [dropdownList, setDropdownList] = useState(props.list ?? []);
  const [isDropdownOpen, setIsDropdownOpen] = useState<boolean>(false);
  const [inputValue, setInputValue] = useState("");

  const dropdownRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    //close dropdown on click outside
    function handleClickOutside(event) {
      if (dropdownRef.current && !dropdownRef.current.contains(event.target)) {
        setIsDropdownOpen(false);
      }
    }

    document.addEventListener("mousedown", handleClickOutside);
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [dropdownRef]);

  useEffect(() => {
    setDropdownList(props.list);
  }, [props.list]);

  useEffect(() => {
    setDropdownList(props.list);
  }, [props.list]);

  console.log(inputValue, "pppp");

  return (
    <div className={styles["dropdown-container"]}>
      <div
        className={styles["input-container"]}
        style={
          { "--border-radius": `${props.borderRadius ?? 8}px` } as CSSProperties
        }
      >
        <input
          disabled={props.isDisabled}
          className={styles.input}
          type="text"
          value={inputValue}
          placeholder={props.inputPlaceholder}
          onFocus={(e) => {
            setIsDropdownOpen(
              e.currentTarget.value.length > 2 || props.openOnClick
            );

            props.onInputFocusCallback &&
              props.onInputFocusCallback(e.currentTarget.value);
          }}
          onChange={(e) => {
            setInputValue(e.currentTarget.value);

            setIsDropdownOpen(
              e.currentTarget.value.length > 2 || props.openOnClick
            );

            props.onInputChangeCallback &&
              props.onInputChangeCallback(e.currentTarget.value);
          }}
          onClick={() => {
            if (props.openOnClick) setIsDropdownOpen(true);
          }}
        />
        <Icon
          name="chev-down"
          fill="#fff"
          className={styles.chev}
          width={32}
          height={32}
          onClick={() => {
            if (props.openOnClick && !props.isDisabled) setIsDropdownOpen(true);
          }}
        />
      </div>
      {isDropdownOpen && (
        <div
          className={`${styles["dropdown-content-container"]} ${props.dropdownContentContainerClassname}`}
          ref={dropdownRef}
        >
          {dropdownList.length === 0 ? (
            <p>No results found.</p>
          ) : (
            dropdownList.map((elt) => {
              return (
                <p
                  onClick={() => {
                    props.onEltClickCallback &&
                      props.onEltClickCallback(elt.elt);

                    setInputValue("");

                    setIsDropdownOpen(false);
                  }}
                >
                  {elt.displayElt}
                </p>
              );
            })
          )}
        </div>
      )}
    </div>
  );
};

Dropdown.defaultProps = {
  openOnClick: false,
};

export default Dropdown;
