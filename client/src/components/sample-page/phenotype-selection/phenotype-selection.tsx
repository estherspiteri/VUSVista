import React, { useState } from "react";
import { IHPOTerm } from "../../../services/sample/sample.dto";
import { SampleService } from "../../../services/sample/sample.service";
import Dropdown from "../../../atoms/dropdown/dropdown";

type PhenotypeSelectionProps = {
  isDisabled?: boolean;
  sampleService?: SampleService;
  onTermClickCallback?: (term: IHPOTerm) => void;
};

const PhenotypeSelection: React.FunctionComponent<PhenotypeSelectionProps> = (
  props: PhenotypeSelectionProps
) => {
  const [HPOTerms, setHPOTerms] = useState<IHPOTerm[]>([]);

  return (
    <Dropdown
      isDisabled={props.isDisabled}
      inputPlaceholder="Type in a phenotype . . ."
      borderRadius={0}
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
      showChev={false}
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
