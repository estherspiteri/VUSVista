import React, { useEffect, useRef, useState } from "react";
import { IHPOTerm } from "../../../services/sample/sample.dto";
import { SampleService } from "../../../services/sample/sample.service";
import Dropdown from "../../../atoms/dropdown/dropdown";

type PhenotypeSelectionProps = {
  sampleService?: SampleService;
  onTermClickCallback?: (term: IHPOTerm) => void;
};

const PhenotypeSelection: React.FunctionComponent<PhenotypeSelectionProps> = (
  props: PhenotypeSelectionProps
) => {
  const [HPOTerms, setHPOTerms] = useState<IHPOTerm[]>([]);

  return (
    <Dropdown
      inputPlaceholder="Type in a phenotype . . ."
      list={HPOTerms.map((t) => {
        return {
          elt: t,
          displayElt: (
            <span>
              {t.ontologyId}: {t.name}
            </span>
          ),
        };
      })}
      onInputFocusCallback={retrieveHPOTerms}
      onInputChangeCallback={retrieveHPOTerms}
      onEltClickCallback={(elt) => {
        props.onTermClickCallback && props.onTermClickCallback(elt as IHPOTerm);
      }}
    />
  );

  function retrieveHPOTerms(term: string) {
    if (term.length > 2) {
      props.sampleService
        ?.getHPOTerms({
          phenotype: term,
        })
        .then((res) => {
          if (res.isSuccess && res.hpoTerms.length > 0) {
            setHPOTerms(res.hpoTerms);
          } else {
            setHPOTerms([]);
          }
        });
    }
  }
};

export default PhenotypeSelection;
