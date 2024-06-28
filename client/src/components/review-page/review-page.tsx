import React, { useState } from "react";
import styles from "./review-page.module.scss";
import {
  ILoadReviewAcmgRules,
  ILoadReviewPublications,
} from "../../models/load-review.model";
import { ReviewService } from "../../services/review/review.service";
import Dropdown from "../../atoms/dropdown/dropdown";
import Icon from "../../atoms/icons/icon";
import TextArea from "../../atoms/text-area/text-area";
import { IVUSSummary } from "../../models/vus-summary.model";
import VariantSummary from "../shared/variant-summary/variant-summary";
import Button from "../../atoms/button/button";
import Modal from "../../atoms/modal/modal";
import Loader from "../../atoms/loader/loader";

type ReviewPageProps = {
  variantSummary: IVUSSummary;
  acmgRules?: ILoadReviewAcmgRules[];
  publications: ILoadReviewPublications[];
  classifications: string[];
  selectedAcmgRuleToAdd?: ILoadReviewAcmgRules;
  selectedAcmgRuleToRemove?: ILoadReviewAcmgRules;
  reviewService?: ReviewService;
};

const ReviewPage: React.FunctionComponent<ReviewPageProps> = (
  props: ReviewPageProps
) => {
  const [selectedClassification, setSelectedClassification] =
    useState<string>("");
  const [typedClassification, setTypedClassification] = useState("");

  const [selectedAcmgRules, setSelectedAcmgRules] = useState<
    ILoadReviewAcmgRules[]
  >(
    props.selectedAcmgRuleToAdd
      ? [props.selectedAcmgRuleToAdd]
      : props.selectedAcmgRuleToRemove
      ? [props.selectedAcmgRuleToRemove]
      : []
  );
  const [typedAcmgRule, setTypedAcmgRule] = useState("");

  const [selectedPublications, setSelectedPublications] = useState<
    ILoadReviewPublications[]
  >([]);
  const [typedPublication, setTypedPublication] = useState("");

  const [reason, setReason] = useState("");

  const [isConfirmSaveModalVisible, setIsConfirmSaveModalVisible] =
    useState(false);

  const [isSavingReview, setIsSavingReview] = useState(false);

  return isSavingReview ? (
    <Loader />
  ) : (
    <div className={styles["review-page-container"]}>
      <div className={styles["title-container"]}>
        <div className={styles.title}>Vus Classification Review</div>
        <div className={styles.description}>
          <p
            dangerouslySetInnerHTML={{
              __html: `Update the variant's classification based on ${
                props.selectedAcmgRuleToAdd || props.selectedAcmgRuleToRemove
                  ? `the <b>${
                      props.selectedAcmgRuleToAdd ? "newly added" : "removed"
                    } ACMG rule ${
                      props.selectedAcmgRuleToAdd
                        ? props.selectedAcmgRuleToAdd.name
                        : props.selectedAcmgRuleToRemove.name
                    }</b>. You may use other `
                  : ""
              } relevant publications ${
                props.selectedAcmgRuleToAdd || props.selectedAcmgRuleToRemove
                  ? `that the variant already has to support this ${
                      props.selectedAcmgRuleToAdd
                        ? "new ACMG rule addition"
                        : "ACMG rule removal"
                    }`
                  : "and/or ACMG rules"
              }. You can also choose to write down ${
                props.selectedAcmgRuleToAdd || props.selectedAcmgRuleToRemove
                  ? "any comments"
                  : "the motivation behind this classification review"
              } in the text box below.`,
            }}
          />
          {!props.selectedAcmgRuleToAdd && !props.selectedAcmgRuleToRemove && (
            <p>
              For a review to be valid at least one of the following must be
              selected or populated: ACMG rule/s, publication/s, written-down
              reason.
            </p>
          )}
        </div>
      </div>
      <div className={styles["variant-info"]}>
        <div className={styles["variant-summary"]}>
          <VariantSummary variant={props.variantSummary} />
        </div>
        <div className={styles["prev-classification-container"]}>
          <p className={styles["section-title"]}>
            This variant's current classification is:
          </p>
          <p
            className={`${styles["prev-classification"]} ${
              styles[
                props.variantSummary.classification
                  .toLowerCase()
                  .replace("_", "-")
              ]
            }`}
          >
            {props.variantSummary.classification.replace("_", " ")}
          </p>
        </div>
      </div>
      <div className={styles["classification-container"]}>
        <p className={styles["section-title"]}>Select its new classification</p>

        {selectedClassification.length === 0 ? (
          <Dropdown
            inputPlaceholder="Select a new classification ..."
            openOnClick={true}
            dropdownContentContainerClassname={
              styles["dropdown-content-container"]
            }
            list={
              props.classifications
                ?.filter(
                  (c) =>
                    !selectedClassification?.includes(c) &&
                    c
                      ?.toLowerCase()
                      .includes(typedClassification?.toLowerCase())
                )
                .map((c) => {
                  return {
                    elt: c,
                    displayElt: <span>{c.replace("_", " ")}</span>,
                  };
                }) ?? []
            }
            onEltClickCallback={(elt) =>
              setTimeout(() => {
                setSelectedClassification((elt as string).replace("_", " "));
              }, 300)
            }
            onInputChangeCallback={(val) => setTypedClassification(val)}
          />
        ) : (
          <div className={styles.selected}>
            <p>{selectedClassification}</p>
            <Icon
              className={styles.close}
              name="close"
              stroke="#008080"
              onClick={() => setSelectedClassification("")}
            />
          </div>
        )}
      </div>

      <div className={styles["dropdown-containers-wrapper"]}>
        <div
          className={`${styles["dropdown-container"]} ${styles["acmg-container"]}`}
        >
          <p className={styles["section-title"]}>ACMG rules selection</p>
          {!props.selectedAcmgRuleToAdd && !props.selectedAcmgRuleToRemove && (
            <Dropdown
              inputPlaceholder="Select from variant's assigned ACMG rules ..."
              openOnClick={true}
              dropdownContentContainerClassname={
                styles["dropdown-content-container"]
              }
              list={
                props.acmgRules
                  ?.filter(
                    (a) =>
                      !selectedAcmgRules?.includes(a) &&
                      a.name
                        ?.toLowerCase()
                        .includes(typedAcmgRule?.toLowerCase())
                  )
                  .map((r) => {
                    return { elt: r, displayElt: <span>{r.name}</span> };
                  }) ?? []
              }
              onEltClickCallback={(elt) =>
                setTimeout(() => {
                  setSelectedAcmgRules(
                    selectedAcmgRules.concat(elt as ILoadReviewAcmgRules)
                  );
                }, 300)
              }
              onInputChangeCallback={(val) => setTypedAcmgRule(val)}
            />
          )}
          <div className={styles["selection-container"]}>
            <p
              dangerouslySetInnerHTML={{
                __html: `This review is based on the
              <b>${
                props.selectedAcmgRuleToAdd
                  ? "addition of the "
                  : props.selectedAcmgRuleToRemove
                  ? "removal of the "
                  : ""
              }</b>following ACMG rule${
                  props.selectedAcmgRuleToAdd || props.selectedAcmgRuleToRemove
                    ? ""
                    : "s"
                }:`,
              }}
            />
            {selectedAcmgRules.length > 0 ? (
              <div className={styles["selected-values-container"]}>
                {selectedAcmgRules.map((r) => (
                  <div
                    className={`${styles.selected} ${styles["selected-acmg"]}`}
                  >
                    <p>{r.name}</p>
                    {!props.selectedAcmgRuleToAdd &&
                      !props.selectedAcmgRuleToRemove && (
                        <Icon
                          className={styles.close}
                          name="close"
                          stroke="#008080"
                          onClick={() =>
                            setSelectedAcmgRules(
                              selectedAcmgRules.filter(
                                (rule) => rule.id !== r.id
                              )
                            )
                          }
                        />
                      )}
                  </div>
                ))}
              </div>
            ) : (
              <p className={styles["no-selection"]}>No selected ACMG rules</p>
            )}
          </div>
        </div>
        <div className={styles["dropdown-container"]}>
          <p className={styles["section-title"]}>Publications selection</p>
          <Dropdown
            inputPlaceholder="Select from variant's linked publications ..."
            openOnClick={true}
            dropdownContentContainerClassname={
              styles["dropdown-content-container"]
            }
            list={
              props.publications
                ?.filter(
                  (p) =>
                    !selectedPublications?.includes(p) &&
                    (p.title
                      ?.toLowerCase()
                      .includes(typedPublication?.toLowerCase()) ||
                      p.doi
                        ?.toLowerCase()
                        ?.includes(typedPublication?.toLowerCase()))
                )
                .map((p) => {
                  return {
                    elt: p,
                    displayElt: (
                      <div>
                        <p>{p.title}</p>
                        <p className={styles.doi}>DOI: {p.doi}</p>
                      </div>
                    ),
                  };
                }) ?? []
            }
            onEltClickCallback={(elt) =>
              setTimeout(() => {
                setSelectedPublications(
                  selectedPublications.concat(elt as ILoadReviewPublications)
                );
              }, 300)
            }
            onInputChangeCallback={(val) => setTypedPublication(val)}
          />
          <div className={styles["selection-container"]}>
            <p>This review is based on the following publications:</p>
            {selectedPublications.length > 0 ? (
              <div className={styles["selected-values-container"]}>
                {selectedPublications.map((p) => (
                  <div className={styles.selected}>
                    <div className={styles["selected-content"]}>
                      <p>{p.title}</p>
                      <p className={styles.doi}>DOI: {p.doi}</p>
                    </div>
                    <Icon
                      className={styles.close}
                      name="close"
                      stroke="#008080"
                      onClick={() =>
                        setSelectedPublications(
                          selectedPublications.filter((pub) => pub.id !== p.id)
                        )
                      }
                    />
                  </div>
                ))}
              </div>
            ) : (
              <p className={styles["no-selection"]}>No selected publications</p>
            )}
          </div>
        </div>
      </div>
      <div className={styles["reason-container"]}>
        <p className={styles["reason-title"]}>Motivation for review</p>
        <TextArea
          placeholder="Type in any comments to justify this classification review ..."
          value={reason}
          className={styles.reason}
          onChange={(e) => setReason(e.currentTarget.value)}
        />
      </div>
      <Button
        text="Save review"
        className={styles.btn}
        disabled={
          selectedClassification.length === 0 ||
          (selectedAcmgRules.length === 0 &&
            selectedPublications.length === 0 &&
            reason.length === 0)
        }
        onClick={() => setIsConfirmSaveModalVisible(true)}
      />

      {isConfirmSaveModalVisible && (
        <Modal
          modalContainerStyle={styles["confirm-save-modal"]}
          isClosable={true}
          onCloseIconClickCallback={() => setIsConfirmSaveModalVisible(false)}
        >
          <div className={styles["confirm-save-modal-content"]}>
            Are you sure you want to save this Classification Review?
            <div className={styles["option-btns"]}>
              <Button text="Yes, save it" onClick={saveReview} />
              <Button
                text="Continue editing"
                onClick={() => setIsConfirmSaveModalVisible(false)}
              />
            </div>
          </div>
        </Modal>
      )}
    </div>
  );

  function saveReview() {
    setIsSavingReview(true);

    props.reviewService
      .saveClassificationReview({
        vusId: props.variantSummary.id,
        newClassification: selectedClassification,
        reason: reason,
        acmgRuleIds: selectedAcmgRules.map((a) => a.id),
        publicationIds: selectedPublications.map((p) => p.id),
        isNewAcmgAdded: props.selectedAcmgRuleToAdd !== undefined,
        isExistingAcmgRemoved: props.selectedAcmgRuleToRemove !== undefined,
      })
      .then((res) => {
        if (res.isSuccess) {
          window.location.href = `/vus/${props.variantSummary.id}`;
        } else {
          setIsSavingReview(false);
        }
      });
  }
};

export default ReviewPage;
