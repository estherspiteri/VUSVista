import React, { useState } from "react";
import styles from "./sample-info.module.scss";
import { Genotype, ISample } from "../../../models/view-samples.model";
import { SampleService } from "../../../services/sample/sample.service";
import VariantSummary from "../../shared/variant-summary/variant-summary";
import Icon from "../../../atoms/icons/icon";
import Text from "../../../atoms/text/text";
import { openInNewWindow } from "../../../helpers/open-links";
import Modal from "../../../atoms/modal/modal";
import Button from "../../../atoms/button/button";
import VusTable from "../../view-all-vus-page/vus-table/vus-table";
import { IVariantToAddInfo } from "../../../models/variant-to-add-info.model";
import Loader from "../../../atoms/loader/loader";
import SamplePhenotypeSelection from "../../shared/sample-phenotype-selection/sample-phenotype-selection";

type SampleInfoProps = {
  sample: ISample;
  sampleService: SampleService;
};

const SampleInfo: React.FunctionComponent<SampleInfoProps> = (
  props: SampleInfoProps
) => {
  const [variants, setVariants] = useState(props.sample.variants);
  const [notSampleVariants, setNotSampleVariants] = useState(
    props.sample.notSampleVariants
  );

  const [variantBeingEdited, setVariantBeingEdited] = useState(undefined);
  const [editedHgvs, setEditedHgvs] = useState(undefined);
  const [isUpdatingHgvs, setIsUpdatingHgvs] = useState(false);
  const [isHgvsUpdateWarningModalVisible, setIsHgvsUpdateWarningModalVisible] =
    useState(false);

  const [isAddingVariantsModalVisible, setIsAddingVariantsModalVisible] =
    useState(false);
  const [variantIdsToAdd, setVariantIdsToAdd] = useState<number[]>([]);
  const [variantInfoToAdd, setVariantInfoToAdd] = useState<IVariantToAddInfo[]>(
    []
  );
  const [showVariantInfoToAdd, setShowVariantInfoToAdd] = useState(false);
  const [isAddingVariants, setIsAddingVariants] = useState(false);

  const [variantIdsToRemove, setVariantIdsToRemove] = useState<number[]>([]);
  const [isRemovingVariantsModalVisible, setIsRemovingVariantsModalVisible] =
    useState(false);
  const [isRemovingVariants, setIsRemovingVariants] = useState(false);
  const [isDeletingSample, setIsDeletingSample] = useState(false);

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
                <div className={styles["delete-variant"]}>
                  <Icon name="bin" fill="#fff" width={21} height={21} />
                </div>
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

                        <input
                          type="checkbox"
                          className={styles.checkbox}
                          checked={variantIdsToRemove.some(
                            (id) => id === v.variantId
                          )}
                          name={v.variantId.toString()}
                          onChange={() => {
                            if (
                              variantIdsToRemove.some(
                                (id) => id === v.variantId
                              )
                            ) {
                              setVariantIdsToRemove(
                                variantIdsToRemove.filter(
                                  (id) => id !== v.variantId
                                )
                              );
                            } else {
                              setVariantIdsToRemove(
                                variantIdsToRemove.concat(v.variantId)
                              );
                            }
                          }}
                        />
                      </div>
                    </div>
                  );
                })}
              </div>
            </div>
          </div>

          {variantIdsToRemove.length > 0 && (
            <Button
              text="Remove selected variant/s"
              icon="bin"
              className={styles["remove-variant-button"]}
              onClick={() =>
                variantIdsToRemove.length === variants.length
                  ? setIsDeletingSample(true)
                  : setIsRemovingVariantsModalVisible(true)
              }
            />
          )}

          {notSampleVariants.length > 0 && (
            <Button
              text="Add variants"
              icon="add"
              className={styles["add-variant-button"]}
              onClick={() => setIsAddingVariantsModalVisible(true)}
            />
          )}
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

      {(isRemovingVariantsModalVisible || isDeletingSample) && (
        <Modal
          title="Please Note"
          isClosable={!isRemovingVariants}
          modalContainerStyle={styles["remove-variants-modal"]}
          onCloseIconClickCallback={closeVariantRemovalModal}
        >
          <div className={styles["remove-variants-modal-content"]}>
            {isDeletingSample ? (
              <p>
                <p>
                  When removing all of the sample's variants, the sample itself
                  is deleted. Are you sure you want to delete sample&nbsp;
                  <b>{props.sample.sampleId}</b> ?
                </p>
              </p>
            ) : (
              <p>
                Are you sure you want to remove the variant/s with the following
                Ids:&nbsp;
                <b>{variantIdsToRemove.join(", ")}</b> from sample&nbsp;
                <b>{props.sample.sampleId}</b> ?
              </p>
            )}
            <div className={styles["option-btns"]}>
              <Button
                disabled={isRemovingVariants}
                text={
                  isDeletingSample
                    ? "Yes, delete sample"
                    : "Yes, remove variant/s"
                }
                onClick={() => removeVariants()}
              />
              <Button
                disabled={isRemovingVariants}
                text="No, return to sample page"
                onClick={closeVariantRemovalModal}
              />
            </div>
          </div>
          {isRemovingVariants && <Loader />}
        </Modal>
      )}

      {isAddingVariantsModalVisible && (
        <Modal
          title={showVariantInfoToAdd ? "Input Variant Info" : "Add Variants"}
          isClosable={!isAddingVariants}
          modalContainerStyle={styles["add-variants-modal"]}
          onCloseIconClickCallback={closeAddVariantsModal}
        >
          <div className={styles["add-variants-modal-content"]}>
            {showVariantInfoToAdd ? (
              <p>
                For each of the selected variants to be added to this sample,
                choose the genotype and input the HGVS.
              </p>
            ) : (
              <p>
                Select which of the below variants you would like to add to
                Sample&nbsp;
                <b>{props.sample.sampleId}</b> by ticking the respective
                checkboxes.
              </p>
            )}

            {showVariantInfoToAdd ? (
              <div className={styles["not-sample-variants"]}>
                {notSampleVariants
                  .filter((v) => variantIdsToAdd.includes(v.variantId))
                  .map((v) => {
                    const variantInfo =
                      variantInfoToAdd.find(
                        (variant) => variant.variantId === v.variantId
                      ) ?? null;

                    return (
                      <div className={styles["variant-to-add"]}>
                        <div className={styles.summary}>
                          <VariantSummary variant={v.variant} />
                        </div>
                        <div className={styles.info}>
                          <span>Genotype:</span>
                          <div className={styles.pills}>
                            {["Heterozygous", "Homozygous"].map((g) => {
                              return (
                                <div
                                  className={`${styles.pill} ${
                                    variantInfo?.genotype === g
                                      ? styles["selected-pill"]
                                      : ""
                                  }`}
                                  onClick={() => {
                                    if (!isAddingVariants) {
                                      updateAddVariantGenotype(
                                        g,
                                        v.variantId,
                                        variantInfo
                                      );
                                    }
                                  }}
                                >
                                  {g}
                                </div>
                              );
                            })}
                          </div>
                        </div>
                        <div className={styles.info}>
                          <span>HGVS:</span>
                          <Text
                            disabled={isAddingVariants}
                            onChange={(e) =>
                              updateAddVariantHgvs(
                                e.currentTarget.value,
                                v.variantId,
                                variantInfo
                              )
                            }
                          />
                        </div>
                      </div>
                    );
                  })}
              </div>
            ) : (
              <VusTable
                vusList={notSampleVariants.map((v) => v.variant)}
                isClickable={false}
                showCheckboxes={true}
                onSelectedVariantsUpdate={(variantsToAdd) =>
                  setVariantIdsToAdd(variantsToAdd)
                }
              />
            )}
            <div className={styles["option-btns"]}>
              {variantIdsToAdd.length > 0 && (
                <Button
                  disabled={
                    (showVariantInfoToAdd &&
                      (variantInfoToAdd.length !== variantIdsToAdd.length ||
                        variantInfoToAdd.filter((v) => !v.hgvs || !v.genotype)
                          .length > 0)) ||
                    isAddingVariants
                  }
                  text={
                    showVariantInfoToAdd
                      ? "Add variants"
                      : "Input variants' info"
                  }
                  onClick={() => {
                    if (variantInfoToAdd.length === 0) {
                      setShowVariantInfoToAdd(true);
                    } else {
                      addVariants();
                    }
                  }}
                />
              )}
              <Button
                text="Return to sample page"
                onClick={closeAddVariantsModal}
                disabled={isAddingVariants}
              />
            </div>
          </div>
          {isAddingVariants && <Loader />}
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

  function closeVariantRemovalModal() {
    setIsDeletingSample(false);
    setIsRemovingVariants(false);
    setIsRemovingVariantsModalVisible(false);
    setVariantIdsToRemove([]);
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

  function closeAddVariantsModal() {
    setIsAddingVariantsModalVisible(false);
    setVariantIdsToAdd([]);
    setShowVariantInfoToAdd(false);
    setVariantInfoToAdd([]);
  }

  function updateAddVariantGenotype(
    genotype: string,
    variantId: number,
    variantInfo?: IVariantToAddInfo
  ) {
    if (variantInfo) {
      let variantsInfoToAddUpdated = [];

      variantInfoToAdd.forEach((info) => {
        if (info.variantId === variantId) {
          variantsInfoToAddUpdated = variantsInfoToAddUpdated.concat({
            ...info,
            genotype: genotype,
          });
        } else {
          variantsInfoToAddUpdated = variantsInfoToAddUpdated.concat(info);
        }
      });

      setVariantInfoToAdd(variantsInfoToAddUpdated);
    } else {
      setVariantInfoToAdd(
        variantInfoToAdd.concat({
          variantId: variantId,
          genotype: genotype,
        })
      );
    }
  }

  function updateAddVariantHgvs(
    hgvs: string,
    variantId: number,
    variantInfo?: IVariantToAddInfo
  ) {
    if (variantInfo) {
      let variantsInfoToAddUpdated = [];

      variantInfoToAdd.forEach((info) => {
        if (info.variantId === variantId) {
          variantsInfoToAddUpdated = variantsInfoToAddUpdated.concat({
            ...info,
            hgvs: hgvs,
          });
        } else {
          variantsInfoToAddUpdated = variantsInfoToAddUpdated.concat(info);
        }
      });

      setVariantInfoToAdd(variantsInfoToAddUpdated);
    } else {
      setVariantInfoToAdd(
        variantInfoToAdd.concat({
          variantId: variantId,
          hgvs: hgvs,
        })
      );
    }
  }

  function addVariants() {
    setIsAddingVariants(true);

    props.sampleService
      .addVariants({
        sampleId: props.sample.sampleId,
        variantsToAdd: variantInfoToAdd,
      })
      .then((res) => {
        if (res.isSuccess) {
          setVariants(res.updatedVariants);
          setNotSampleVariants(res.updatedNotSampleVariants);
          closeAddVariantsModal();
        }
        setIsAddingVariants(false);
      });
  }

  function removeVariants() {
    setIsRemovingVariants(true);

    props.sampleService
      .removeVariants({
        sampleId: props.sample.sampleId,
        variantIdsToRemove: variantIdsToRemove,
        isDeleteSample: isDeletingSample,
      })
      .then((res) => {
        if (res.isSuccess) {
          if (res.isSampleDeleted) {
            window.location.href = "/view-samples";
          } else {
            setVariants(res.updatedVariants);
            setNotSampleVariants(res.updatedNotSampleVariants);
            closeVariantRemovalModal();
          }
        }
        setIsRemovingVariants(false);
      });
  }
};

export default SampleInfo;
