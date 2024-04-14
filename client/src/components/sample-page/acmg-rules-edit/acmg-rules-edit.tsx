import React, { useState } from "react";
import styles from "./acmg-rules-edit.module.scss";
import Icon from "../../../atoms/icons/icon";
import { IAcmgRule } from "../../../models/acmg-rule.model";
import { VusService } from "../../../services/vus/vus.service";

type AcmgRulesEditProps = {
  variantId?: number;
  variantAcmgRuleIds?: number[];
  allAcmgRules: IAcmgRule[];
  isAcmgMenuClosable?: boolean;
  vusService?: VusService;
  onEditIconClick?: (isOpen: boolean) => void;
  onMenuAcmgRuleHover?: (amcgRuleId?: number) => void;
  onAcmgRulesSelectionUpdate?: (acmgRuleIds?: number[]) => void;
};

const AcmgRulesEdit: React.FunctionComponent<AcmgRulesEditProps> = (
  props: AcmgRulesEditProps
) => {
  const [isAddMenuVisible, setIsAddMenuVisible] = useState(
    !props.isAcmgMenuClosable
  );
  const [selectedAcmgRules, setSelectedAcmgRules] = useState(
    props.variantAcmgRuleIds ?? []
  );
  const [mouseOverRule, setMouseOverRule] = useState<number>(undefined);

  return (
    <div className={styles.acmg}>
      <div className={styles["selected-rules"]}>
        {/** Selected Acmg Rules */}
        {props.allAcmgRules
          .filter((r) => selectedAcmgRules.includes(r.id))
          ?.map((r) => {
            const styleClass = r.name
              .substring(0, r.name.length - 1)
              .toLowerCase();

            return (
              <div
                onMouseOver={() => setMouseOverRule(r.id)}
                onMouseLeave={() => setMouseOverRule(undefined)}
                className={`${styles["selected-acmg"]} ${
                  styles[`${styleClass}`]
                }`}
              >
                {mouseOverRule === r.id ? (
                  <Icon
                    name="bin"
                    className={styles.bin}
                    onClick={() => removeAcmgRule(r.id)}
                  />
                ) : (
                  r.name
                )}
              </div>
            );
          })}
      </div>

      {/** Add Menu */}
      {getAvailableAcmgRules().length > 0 &&
        (isAddMenuVisible ? (
          <div className={styles["add-menu"]}>
            {props.isAcmgMenuClosable && (
              <Icon
                name="remove"
                fill="white"
                width={20}
                height={20}
                className={styles["add-menu-icon"]}
                onClick={() => {
                  setIsAddMenuVisible(false);
                  props.onEditIconClick && props.onEditIconClick(false);
                }}
              />
            )}
            {getAvailableAcmgRules().map((r) => {
              return (
                <div
                  onClick={() => addAcmgRule(r.id)}
                  onMouseOver={() => {
                    props.onMenuAcmgRuleHover &&
                      props.onMenuAcmgRuleHover(r.id);
                  }}
                  onMouseLeave={() => {
                    props.onMenuAcmgRuleHover &&
                      props.onMenuAcmgRuleHover(undefined);
                  }}
                >
                  {r.name}
                </div>
              );
            })}
          </div>
        ) : (
          <Icon
            className={styles["add-acmg"]}
            name="add-outline"
            fill="white"
            width={35}
            height={35}
            onClick={() => {
              setIsAddMenuVisible(true);
              props.onEditIconClick && props.onEditIconClick(true);
            }}
          />
        ))}
    </div>
  );

  function addAcmgRule(rule_id: number) {
    const updatedAcmgRules = selectedAcmgRules.concat(rule_id);
    setSelectedAcmgRules(updatedAcmgRules);

    if (props.variantId) {
      props.vusService?.addAcmgRule({
        variantId: props.variantId,
        ruleId: rule_id,
      });
    }

    props.onAcmgRulesSelectionUpdate &&
      props.onAcmgRulesSelectionUpdate(updatedAcmgRules);
  }

  function removeAcmgRule(rule_id: number) {
    const updatedAcmgRules = selectedAcmgRules.filter((id) => id !== rule_id);
    setSelectedAcmgRules(updatedAcmgRules);

    if (props.variantId) {
      props.vusService?.removeAcmgRule({
        variantId: props.variantId,
        ruleId: rule_id,
      });
    }

    props.onAcmgRulesSelectionUpdate &&
      props.onAcmgRulesSelectionUpdate(updatedAcmgRules);
  }

  function getAvailableAcmgRules() {
    const availableAcmgRules = props.allAcmgRules.filter(
      (r) => !selectedAcmgRules?.includes(r.id)
    );

    return availableAcmgRules;
  }
};

AcmgRulesEdit.defaultProps = {
  isAcmgMenuClosable: true,
};

export default AcmgRulesEdit;
