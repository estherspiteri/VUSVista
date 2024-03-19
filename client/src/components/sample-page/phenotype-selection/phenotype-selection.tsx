import React, { useEffect, useRef, useState } from "react";
import styles from "./phenotype-selection.module.scss";
import { IHPOTerm } from "../../../services/sample/sample.dto";
import { SampleService } from "../../../services/sample/sample.service";

type PhenotypeSelectionProps = {
  sampleService?: SampleService;
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
        placeholder="Type in a phenotype . . ."
        onFocus={retrieveHPOTerms}
        onChange={retrieveHPOTerms}
      />
      {showHPOTerms && (
        <div className={styles["dropdown-container"]} ref={dropdownRef}>
          {HPOTerms.length === 0 ? (
            <p>No results found.</p>
          ) : (
            HPOTerms.map((term) => {
              return (
                <p
                  onClick={() => {
                    props.onTermClickCallback &&
                      props.onTermClickCallback(term);
                    setShowHPOTerms(false);
                  }}
                >
                  {term.ontologyId}: {term.name}
                </p>
              );
            })
          )}
        </div>
      )}
    </div>
  );

  function retrieveHPOTerms(e) {
    if (e.currentTarget.value.length > 2) {
      props.sampleService
        ?.getHPOTerms({
          phenotype: e.currentTarget?.value,
        })
        .then((res) => {
          if (res.isSuccess && res.hpoTerms.length > 0) {
            setHPOTerms(res.hpoTerms);
          } else {
            setHPOTerms([]);
          }
          setShowHPOTerms(true);
        });
    }
  }
};

export default PhenotypeSelection;
