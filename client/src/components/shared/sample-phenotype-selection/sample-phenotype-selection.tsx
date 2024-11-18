import React, { useState } from "react";
import styles from "./sample-phenotype-selection.module.scss";
import Icon from "../../../atoms/icons/icon";
import { IHPOTerm } from "../../../services/sample/sample.dto";
import { SampleService } from "../../../services/sample/sample.service";
import { openInNewWindow } from "../../../helpers/open-links";
import PhenotypeSelection from "../../sample-page/phenotype-selection/phenotype-selection";

type SamplePhenotypeSelectionProps = {
  sampleId?: string;
  isDisabled?: boolean;
  selectedPhenotypes?: IHPOTerm[];
  sampleService?: SampleService;
  isSelectingPhenotypesForVariant?: boolean;
  onPhenotypesUpdateCallback?: (phenotypes: IHPOTerm[]) => void;
};

const SamplePhenotypeSelection: React.FunctionComponent<
  SamplePhenotypeSelectionProps
> = (props: SamplePhenotypeSelectionProps) => {
  const [samplePhenotypesSelection, setSamplePhenotypesSelection] = useState<
    IHPOTerm[]
  >(props.selectedPhenotypes ?? []);

  return (
    <div className={styles["sample-phenotype-selection-container"]}>
      <PhenotypeSelection
        sampleService={props.sampleService}
        onTermClickCallback={(term) => addToSelection(term)}
      />
      <div className={styles["phenotypes-selected"]}>
        {samplePhenotypesSelection !== undefined &&
          samplePhenotypesSelection.length > 0 &&
          samplePhenotypesSelection.map((term) => {
            return (
              <p className={styles.phenotype}>
                <span
                  className={styles["phenotype-term"]}
                  onClick={() =>
                    !props.isDisabled &&
                    openInNewWindow(
                      `https://hpo.jax.org/app/browse/term/${term.ontologyId}`
                    )
                  }
                >
                  <b>{term.ontologyId}</b>: {term.name}
                </span>
                <Icon
                  name="close"
                  onClick={() => !props.isDisabled && removeSelection(term)}
                  stroke="#008080"
                />
              </p>
            );
          })}
      </div>
    </div>
  );

  function removeSelection(term: IHPOTerm) {
    if (samplePhenotypesSelection.some((s) => s === term)) {
      const updatedSelection = samplePhenotypesSelection.filter(
        (s) => s !== term
      );

      if (!props.isSelectingPhenotypesForVariant) {
        props.sampleService.removePhenotype({
          sampleId: props.sampleId,
          phenotype: term,
        });
      } else if (props.onPhenotypesUpdateCallback) {
        props.onPhenotypesUpdateCallback(updatedSelection);
      }

      setSamplePhenotypesSelection(updatedSelection);
    }
  }

  function addToSelection(term: IHPOTerm) {
    var updatedSelection = samplePhenotypesSelection;

    const isHPOTermSelectedAlready = samplePhenotypesSelection
      .map((s) => s.ontologyId)
      .includes(term.ontologyId);

    if (!isHPOTermSelectedAlready) {
      updatedSelection = samplePhenotypesSelection.concat([term]);

      if (!props.isSelectingPhenotypesForVariant) {
        // add phenotype to db
        props.sampleService.addPhenotype({
          sampleId: props.sampleId,
          phenotype: term,
        });
      } else if (props.onPhenotypesUpdateCallback) {
        props.onPhenotypesUpdateCallback(updatedSelection);
      }
    }

    setSamplePhenotypesSelection(updatedSelection);
  }
};

SamplePhenotypeSelection.defaultProps = {
  isSelectingPhenotypesForVariant: false,
  isDisabled: false,
};

export default SamplePhenotypeSelection;
