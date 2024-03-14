import React, { useState } from "react";
import styles from "./acmg-rules-edit.module.scss";
import Icon from "../../../atoms/icon/icon";
import { SampleService } from "../../../services/sample/sample.service";

type AcmgRulesEditProps = {
  sampleId: string;
  variantId: number;
  variantAcmgRules: string[];
  allAcmgRules: string[];
  sampleService: SampleService;
};

const AcmgRulesEdit: React.FunctionComponent<AcmgRulesEditProps> = (
  props: AcmgRulesEditProps
) => {
  const [isAddMenuVisible, setIsAddMenuVisible] = useState(false);
  const [selectedAcmgRules, setSelectedAcmgRules] = useState(
    props.variantAcmgRules ?? []
  );
  const [mouseOverRule, setMouseOverRule] = useState("");

  return (
    <div className={styles.acmg}>
      {/** Selected Acmg Rules */}
      {selectedAcmgRules?.map((r) => {
        const styleClass = r.substring(0, r.length - 1).toLowerCase();

        return (
          <div
            onMouseOver={() => setMouseOverRule(r)}
            onMouseLeave={() => setMouseOverRule("")}
            className={`${styles["selected-acmg"]} ${styles[`${styleClass}`]}`}
          >
            {mouseOverRule === r ? (
              <Icon
                name="bin"
                className={styles.bin}
                onClick={() => removeAcmgRule(r)}
              />
            ) : (
              r
            )}
          </div>
        );
      })}

      {/** Add Menu */}
      {isAddMenuVisible ? (
        <div className={styles["add-menu"]}>
          <Icon
            name="add"
            fill="white"
            width={20}
            height={20}
            className={styles["add-menu-icon"]}
            onClick={() => setIsAddMenuVisible(false)}
          />
          {getAvailableAcmgRules().map((r) => {
            return <div onClick={() => addAcmgRule(r)}>{r}</div>;
          })}
        </div>
      ) : (
        <Icon
          className={styles["add-acmg"]}
          name="add-outline"
          fill="white"
          width={37}
          height={37}
          onClick={() => {
            if (getAvailableAcmgRules().length > 0) {
              setIsAddMenuVisible(true);
            }
          }}
        />
      )}
    </div>
  );

  function addAcmgRule(rule_name: string) {
    setSelectedAcmgRules(selectedAcmgRules.concat(rule_name));

    props.sampleService.addAcmgRule({
      sampleId: props.sampleId,
      variantId: props.variantId,
      ruleName: rule_name,
    });
  }

  function removeAcmgRule(rule_name: string) {
    setSelectedAcmgRules(selectedAcmgRules.filter((r) => r !== rule_name));

    props.sampleService.removeAcmgRule({
      sampleId: props.sampleId,
      variantId: props.variantId,
      ruleName: rule_name,
    });
  }

  function getAvailableAcmgRules() {
    const availableAcmgRules = props.allAcmgRules.filter(
      (r) => !selectedAcmgRules?.includes(r)
    );

    return availableAcmgRules;
  }
};

export default AcmgRulesEdit;
