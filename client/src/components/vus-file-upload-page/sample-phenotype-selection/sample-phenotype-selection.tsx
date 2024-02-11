import React, { useState } from "react";
import styles from "./sample-phenotype-selection.module.scss";
import { VusService } from "../../../services/vus/vus.service";
import Icon from "../../../atoms/icon/icon";
import { IHPOTerm } from "../../../services/vus/vus.dto";
import PhenotypeSelection from "../phenotype-selection/phenotype-selection";

type SamplePhenotypeSelectionProps = {
  sampleId: string;
  vusService?: VusService;
  onSamplePhenotypesSelectionUpdate?: (
    samplePhenotypesSelection: IHPOTerm[]
  ) => void;
};

const SamplePhenotypeSelection: React.FunctionComponent<
  SamplePhenotypeSelectionProps
> = (props: SamplePhenotypeSelectionProps) => {
  const [samplePhenotypesSelection, setSamplePhenotypesSelection] = useState<
    IHPOTerm[]
  >([]);

  return (
    <div className={styles["sample-phenotype-selection-container"]}>
      <div className={styles["sample-id-container"]}>
        <div>{props.sampleId}</div>
      </div>
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
                <span>
                  <b>{term.ontologyId}</b>: {term.name}
                </span>
              </p>
            );
          })}
      </div>
      <div className={styles['phenotype-selection']}>
        <PhenotypeSelection
          vusService={props.vusService}
          onTermClickCallback={(term) => addToSelection(term)}
        />
      </div>
    </div>
  );

  function removeSelection(term: IHPOTerm) {
    const updatedSelection = samplePhenotypesSelection.filter(
      (s) => s !== term
    );

    setSamplePhenotypesSelection(updatedSelection);

    props.onSamplePhenotypesSelectionUpdate &&
      props.onSamplePhenotypesSelectionUpdate(updatedSelection);
  }

  function addToSelection(term: IHPOTerm) {
    var updatedSelection = samplePhenotypesSelection;

    const isHPOTermSelectedAlready = samplePhenotypesSelection
      .map((s) => s.ontologyId)
      .includes(term.ontologyId);

    if (!isHPOTermSelectedAlready) {
      updatedSelection = samplePhenotypesSelection.concat([term]);
    }

    setSamplePhenotypesSelection(updatedSelection);

    props.onSamplePhenotypesSelectionUpdate &&
      props.onSamplePhenotypesSelectionUpdate(updatedSelection);
  }
};

export default SamplePhenotypeSelection;
