import React from "react";
import styles from "./acmg-rule-info.module.scss";
import { IAcmgRule } from "../../../models/acmg-rule.model";

type AcmgRuleInfoProps = {
  acmgRule?: IAcmgRule;
};

const AcmgRuleInfo: React.FunctionComponent<AcmgRuleInfoProps> = (
  props: AcmgRuleInfoProps
) => {
  return (
    <div className={styles["acmg-rule-info-container"]}>
      {props.acmgRule ? (
        <>
          <p className={styles.title}>{props.acmgRule.name}</p>
          <p className={styles.description}>{props.acmgRule.description}</p>
          <p className={styles.strength}>
            <p>
              Pathogenicity:{" "}
              <b>
                {props.acmgRule.name.startsWith("P") ? "PATHOGENIC" : "BENIGN"}
              </b>
            </p>
            <p>
              Strength: <b>{props.acmgRule.defaultStrength}</b>
            </p>
          </p>
        </>
      ) : (
        <>
          <p className={styles.title}>ACMG Rule Info</p>
          <p className={styles.description}>
            Hover on a rule from the ACMG menu (+) to view information related to the rule.
          </p>
        </>
      )}
    </div>
  );
};

export default AcmgRuleInfo;
