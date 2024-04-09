import React, { useEffect, useRef, useState } from "react";
import styles from "./dropdown.module.scss";

type DropdownElt = { elt: {}; displayElt: JSX.Element };

type DropdownProps = {
  inputPlaceholder?: string;
  list?: { elt: {}; displayElt: JSX.Element }[];
  onEltClickCallback?: (elt: {}) => void;
  onInputFocusCallback?: (val: string) => void;
  onInputChangeCallback?: (val: string) => void;
};

const Dropdown: React.FunctionComponent<DropdownProps> = (
  props: DropdownProps
) => {
  const [dropdownList, setDropdownList] = useState(props.list ?? []);
  const [isDropdownOpen, setIsDropdownOpen] = useState<boolean>(false);

  const dropdownRef = useRef<HTMLDivElement>(null);

  console.log(dropdownList);
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

  return (
    <div className={styles["dropdown-container"]}>
      <input
        className={styles.input}
        type="text"
        placeholder={props.inputPlaceholder}
        onFocus={(e) => {
          setIsDropdownOpen(e.currentTarget.value.length > 2);
          
          props.onInputFocusCallback &&
            props.onInputFocusCallback(e.currentTarget.value);
        }}
        onChange={(e) => {
          setIsDropdownOpen(e.currentTarget.value.length > 2);

          props.onInputChangeCallback &&
            props.onInputChangeCallback(e.currentTarget.value);
        }}
      />
      {isDropdownOpen && (
        <div className={styles["dropdown-container"]} ref={dropdownRef}>
          {dropdownList.length === 0 ? (
            <p>No results found.</p>
          ) : (
            dropdownList.map((elt) => {
              return (
                <p
                  onClick={() => {
                    props.onEltClickCallback &&
                      props.onEltClickCallback(elt.elt);
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

export default Dropdown;
