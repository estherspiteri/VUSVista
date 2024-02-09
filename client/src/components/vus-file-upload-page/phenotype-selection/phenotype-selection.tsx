import React, { useEffect, useRef, useState } from "react";
import { VusService } from "../../../services/vus/vus.service";
import { IHPOTerm } from "../../../services/vus/vus.dto";
import styles from "./phenotype-selection.module.scss";

type PhenotypeSelectionProps = {
  vusService?: VusService;
  onTermClickCallback?: (term: IHPOTerm) => void;
};

const PhenotypeSelection: React.FunctionComponent<PhenotypeSelectionProps> = (
  props: PhenotypeSelectionProps
) => {
  const [HPOTerms, setHPOTerms] = useState<IHPOTerm[]>([]);
  const [showHPOTerms, setShowHPOTerms] = useState<boolean>(false);

  const dropdownRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    //close hpo term options on click outside
    function handleClickOutside(event) {
      if (dropdownRef.current && !dropdownRef.current.contains(event.target)) {
        setShowHPOTerms(false);
      }
    }

    document.addEventListener("mousedown", handleClickOutside);
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [dropdownRef]);

  return (
    <div className={styles["phenotype-selection-input-container"]}>
      <input
        className={styles["phenotype-selection-input"]}
        type="text"
        placeholder="Select phenotype"
        onFocus={retrieveHPOTerms}
        onChange={retrieveHPOTerms}
      />
      {showHPOTerms && (
        <div className={styles["dropdown-container"]} ref={dropdownRef}>
          {HPOTerms.map((term) => {
            return (
              <p
                onClick={() => {
                  props.onTermClickCallback && props.onTermClickCallback(term);
                }}
              >
                {term.ontologyId}: {term.name}
              </p>
            );
          })}
        </div>
      )}
    </div>
  );

  function retrieveHPOTerms(e) {
    if (e.currentTarget.value.length > 3) {
      props.vusService
        ?.getHPOTerms({
          phenotype: e.currentTarget?.value,
        })
        .then((res) => {
          if (res.isSuccess && res.hpoTerms.length > 0) {
            setHPOTerms(res.hpoTerms);
            setShowHPOTerms(true);
          }
        });
    }
  }
};

export default PhenotypeSelection;
