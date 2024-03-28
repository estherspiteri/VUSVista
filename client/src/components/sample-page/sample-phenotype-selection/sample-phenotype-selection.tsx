import React, { useState } from "react";
import styles from "./sample-phenotype-selection.module.scss";
import Icon from "../../../atoms/icons/icon";
import PhenotypeSelection from "../phenotype-selection/phenotype-selection";
import { IHPOTerm } from "../../../services/sample/sample.dto";
import { SampleService } from "../../../services/sample/sample.service";
import { openInNewWindow } from "../../../helpers/open-links";

type SamplePhenotypeSelectionProps = {
  sampleId: string;
  selectedPhenotypes?: IHPOTerm[];
  sampleService?: SampleService;
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
              <p className={styles["phenotype"]}>
                <Icon
                  name="close"
                  onClick={() => removeSelection(term)}
                  stroke="#008080"
                />
                <span
                  className={styles["phenotype-term"]}
                  onClick={() =>
                    openInNewWindow(
                      `https://hpo.jax.org/app/browse/term/${term.ontologyId}`
                    )
                  }
                >
                  <b>{term.ontologyId}</b>: {term.name}
                </span>
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

      props.sampleService.removePhenotype({
        sampleId: props.sampleId,
        phenotype: term,
      });

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

      // add phenotype to db
      props.sampleService.addPhenotype({
        sampleId: props.sampleId,
        phenotype: term,
      });
    }

    setSamplePhenotypesSelection(updatedSelection);
  }
};

export default SamplePhenotypeSelection;
