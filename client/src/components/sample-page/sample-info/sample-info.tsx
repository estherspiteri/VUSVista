import React, { useState } from "react";
import styles from "./sample-info.module.scss";
import { Genotype, ISample } from "../../../models/view-samples.model";
import { SampleService } from "../../../services/sample/sample.service";
import VariantSummary from "../../shared/variant-summary/variant-summary";
import SamplePhenotypeSelection from "../sample-phenotype-selection/sample-phenotype-selection";
import Icon from "../../../atoms/icons/icon";
import Text from "../../../atoms/text/text";
import { openInNewWindow } from "../../../helpers/open-links";
import Modal from "../../../atoms/modal/modal";
import Button from "../../../atoms/button/button";

type SampleInfoProps = {
  sample: ISample;
  sampleService: SampleService;
};

const SampleInfo: React.FunctionComponent<SampleInfoProps> = (
  props: SampleInfoProps
) => {
  const [variants, setVariants] = useState(props.sample.variants);
  const [variantBeingEdited, setVariantBeingEdited] = useState(undefined);
  const [editedHgvs, setEditedHgvs] = useState(undefined);
  const [isUpdatingHgvs, setIsUpdatingHgvs] = useState(false);
  const [isHgvsUpdateWarningModalVisible, setIsHgvsUpdateWarningModalVisible] =
    useState(false);

  return (
    <div className={styles["sample-info-container"]}>
      <div className={styles["sample-info"]}>
        <div className={styles["top-container"]}>
          {/** General information */}
          <div className={styles.information}>
            <div className={styles.info}>
              <div className={styles["info-title"]}>Sample Id:</div>
              {props.sample.sampleId}
            </div>

            <div className={styles.info}>
              <div className={styles["info-title"]}>Genome Version:</div>
              {props.sample.genomeVersion}
            </div>
          </div>
        </div>

        {/** Phenotypes */}
        <div className={styles.phenotypes}>
          <div className={styles["info-title"]}>Phenotypes:</div>
          <SamplePhenotypeSelection
            sampleService={props.sampleService}
            sampleId={props.sample.sampleId}
            selectedPhenotypes={props.sample.phenotype}
          />
        </div>

        {/** Variants */}
        <div className={styles["sample-variant-container"]}>
          <div className={styles["info-title"]}>Variants:</div>
          <div className={styles["sample-variant-container-content"]}>
            <div className={styles["variants-container"]}>
              <div className={styles["variant-titles"]}>
                <div className={styles["variant-summary"]}>Variant Summary</div>
                <div className={styles.genotype}>Genotype</div>
                <div>HGVS</div>
              </div>
              <div className={styles.variants}>
                {variants.map((v) => {
                  return (
                    <div className={styles.variant}>
                      <Icon
                        name="external-link"
                        className={styles["external-link"]}
                        onClick={() => openInNewWindow(`/vus/${v.variantId}`)}
                      />
                      <div className={styles["variant-info-container"]}>
                        <div className={styles["variant-info"]}>
                          <div className={styles["variant-summary"]}>
                            {/*TODO: on click take to VUS page*/}
                            <VariantSummary variant={v.variant} />
                          </div>
                          <div className={styles.genotype}>
                            {v.genotype === Genotype.Heterozygous ? "Aa" : "AA"}
                          </div>
                          <div className={styles.hgvs}>
                            {variantBeingEdited === v.variantId ? (
                              <Text
                                autoFocus={true}
                                value={editedHgvs}
                                disabled={isUpdatingHgvs}
                                onChange={(e) =>
                                  setEditedHgvs(e.currentTarget.value)
                                }
                              />
                            ) : (
                              <span>{v.hgvs}</span>
                            )}
                            <Icon
                              className={styles["hgvs-edit"]}
                              name={
                                variantBeingEdited === v.variantId
                                  ? "save"
                                  : "edit"
                              }
                              fill="#fff"
                              width={16}
                              height={16}
                              onClick={() => {
                                if (!isUpdatingHgvs) {
                                  variantBeingEdited === v.variantId
                                    ? updatedVariantHgvsCheck()
                                    : editVariantHgvs(v.variantId, v.hgvs);
                                }
                              }}
                            />
                            {variantBeingEdited === v.variantId && (
                              <Icon
                                className={styles["hgvs-edit-close"]}
                                name={"close"}
                                stroke="#fff"
                                width={24}
                                height={24}
                                onClick={() => {
                                  setEditedHgvs(undefined);
                                  setVariantBeingEdited(undefined);
                                }}
                              />
                            )}
                            {v.isHgvsUpdated && (
                              <span className={styles["hgvs-updated"]}>
                                Updated
                              </span>
                            )}
                          </div>
                        </div>
                      </div>
                    </div>
                  );
                })}
              </div>
            </div>
          </div>
        </div>
      </div>

      {isHgvsUpdateWarningModalVisible && (
        <Modal
          title="Please Note"
          isClosable={true}
          modalContainerStyle={styles["confirm-update-modal"]}
          onCloseIconClickCallback={closeHgvsCheckModal}
        >
          <div className={styles["confirm-update-modal-content"]}>
            <p>
              Updating an HGVS for this sample's variant will update the HGVS of
              any other samples' variants which have the same HGVS.
            </p>
            <p>Are you sure you want to proceed with the HGVS update?</p>
            <div className={styles["option-btns"]}>
              <Button
                text="Yes, proceed"
                onClick={() =>
                  updateVariantHgvs(variantBeingEdited, editedHgvs)
                }
              />
              <Button
                text="No, return to sample page"
                onClick={closeHgvsCheckModal}
              />
            </div>
          </div>
        </Modal>
      )}
    </div>
  );

  function editVariantHgvs(variantId: number, hgvs: string) {
    setVariantBeingEdited(variantId);
    setEditedHgvs(hgvs);
  }

  function closeHgvsCheckModal() {
    setIsHgvsUpdateWarningModalVisible(false);
    setEditedHgvs(undefined);
    setVariantBeingEdited(undefined);
  }

  function updatedVariantHgvsCheck() {
    //check that hgvs has been updated & make sure it is not an empty string
    if (
      editedHgvs?.trim().length > 0 &&
      variants.find((v) => v.variantId === variantBeingEdited).hgvs !==
        editedHgvs
    ) {
      setIsHgvsUpdateWarningModalVisible(true);
    } else {
      setEditedHgvs(undefined);
      setVariantBeingEdited(undefined);
    }
  }

  function updateVariantHgvs(variantId: number, editedHgvs: any) {
    setIsHgvsUpdateWarningModalVisible(false);
    setIsUpdatingHgvs(true);

    let updatedVariants = [];

    variants.forEach((v) => {
      let variant = v;

      if (v.variantId === variantId) {
        variant.hgvs = editedHgvs;
        variant.isHgvsUpdated = true;
      }

      updatedVariants = updatedVariants.concat(variant);
    });

    props.sampleService
      .updateHgvs({
        variantId: variantId.toString(),
        sampleId: props.sample.sampleId,
        updatedHgvs: editedHgvs,
      })
      .then((res) => {
        if (res.isSuccess) {
          setVariants(updatedVariants);
          setIsUpdatingHgvs(false);
        }
      });

    setEditedHgvs(undefined);
    setVariantBeingEdited(undefined);
  }
};

export default SampleInfo;
