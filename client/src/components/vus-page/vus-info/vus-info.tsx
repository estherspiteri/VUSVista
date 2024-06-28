import React, { useState } from "react";
import styles from "./vus-info.module.scss";
import { IVus } from "../../../models/view-vus.model";
import VariantSummary from "../../shared/variant-summary/variant-summary";
import Icon from "../../../atoms/icons/icon";
import { openInNewWindow } from "../../../helpers/open-links";
import { IAcmgRule } from "../../../models/acmg-rule.model";
import AcmgRulesEdit from "../../sample-page/acmg-rules-edit/acmg-rules-edit";
import { VusService } from "../../../services/vus/vus.service";
import AcmgRuleInfo from "../../sample-page/acmg-rule-info/acmg-rule-info";
import { IClinvarUpdate } from "../../../models/clinvar-updates.model";
import Modal from "../../../atoms/modal/modal";
import Loader from "../../../atoms/loader/loader";
import CalendarDisplay from "../../../atoms/calendar-display/calendar-display";
import Button from "../../../atoms/button/button";
import Text from "../../../atoms/text/text";
import { ISampleToAddInfo } from "../../../models/sample-to-add-info.model";
import SampleTable from "../../view-samples-page/sample-table/sample-table";
import { SampleService } from "../../../services/sample/sample.service";
import SamplePhenotypeSelection from "../../shared/sample-phenotype-selection/sample-phenotype-selection";
import { IPhenotype } from "../../../models/phenotype.model";

type VusInfoProps = {
  vus: IVus;
  acmgRules: IAcmgRule[];
  vusService?: VusService;
  sampleService?: SampleService;
};

const VusInfo: React.FunctionComponent<VusInfoProps> = (
  props: VusInfoProps
) => {
  const [acmgRuleHover, setAcmgRuleHover] = useState<number | undefined>(
    undefined
  );
  const [isAcmgEditMenuOpen, setIsAcmgEditMenuOpen] = useState(false);

  const [showClinvarArchiveModal, setShowClinvarArchiveModal] = useState(false);
  const [clinvarUpdates, setClinvarUpdates] = useState<IClinvarUpdate[]>([]);
  const [datesWithUpdates, setDatesWithUpdates] = useState([]);

  const [phenotypes, setPhenotypes] = useState(props.vus.phenotypes);
  const [samples, setSamples] = useState(props.vus.samples);
  const [notVariantSamples, setNotVariantSamples] = useState(
    props.vus.notVusSamples
  );

  const [isAddingNewSampleModalVisible, setIsAddingNewSampleModalVisible] =
    useState(false);
  const [showNewSampleIdError, setShowNewSampleIdError] = useState(false);

  const [isAddingSamplesModalVisible, setIsAddingSamplesModalVisible] =
    useState(false);
  const [sampleIdsToAdd, setSampleIdsToAdd] = useState<string[]>([]);
  const [sampleInfoToAdd, setSampleInfoToAdd] = useState<ISampleToAddInfo[]>(
    []
  );
  const [showSamplesInfoToAdd, setShowSamplesInfoToAdd] = useState(false);
  const [isAddingSamples, setIsAddingSamples] = useState(false);

  const [sampleIdsToRemove, setSampleIdsToRemove] = useState<string[]>([]);
  const [isRemovingSamplesModalVisible, setIsRemovingSamplesModalVisible] =
    useState(false);
  const [isRemovingSamples, setIsRemovingSamples] = useState(false);

  return (
    <div className={styles["vus-info-container"]}>
      <div className={styles["vus-info"]}>
        <div className={styles["vus-summary"]}>
          <VariantSummary
            variant={{
              id: props.vus.id,
              chromosome: props.vus.chromosome,
              chromosomePosition: props.vus.chromosomePosition,
              gene: props.vus.gene,
              refAllele: props.vus.refAllele,
              altAllele: props.vus.altAllele,
            }}
          />
          <div
            className={`${styles.classification} ${
              styles[props.vus.classification.toLowerCase().replace("_", "-")]
            }`}
          >
            {props.vus.classification.replace("_", " ")}
          </div>
        </div>
        <div className={styles["top-container"]}>
          {/** General information */}
          <div className={styles.information}>
            <div className={styles.info}>
              <div className={styles["info-title"]}>Type:</div>
              {props.vus.type}
            </div>

            <div className={styles.info}>
              <div className={styles["info-title"]}>Homozygotes:</div>
              {props.vus.numHomozygous}
            </div>

            <div className={styles.info}>
              <div className={styles["info-title"]}>Heterozygotes:</div>
              {props.vus.numHeterozygous}
            </div>
          </div>
        </div>

        {/*TODO: Next to is RSID verified do info icon - on hover show what info was compared. Same for clinvar*/}
        {/** External References */}
        <div className={styles["external-ref"]}>
          {/** Clinvar */}
          <div className={styles["info-container"]}>
            <div className={styles["external-ref-title-container"]}>
              <p className={styles["info-title"]}>
                {`ClinVar${
                  props.vus.clinvarClassification?.length > 0 &&
                  !props.vus.rsidDbsnpVerified
                    ? " of suggested dbSNP RSID"
                    : ""
                }:`}
              </p>
              {(props.vus.clinvarClassification?.length > 0 ||
                props.vus.clinvarErrorMsg?.length > 0) && (
                <Icon
                  name="external-link"
                  className={`${styles["external-link"]} ${styles.clinvar} ${
                    props.vus.clinvarClassification?.length === 0
                      ? styles.disabled
                      : props.vus.rsid?.length > 0 &&
                        !props.vus.rsidDbsnpVerified
                      ? styles["unverified-rsid"]
                      : styles.active
                  }`}
                  onClick={(e) => {
                    if (props.vus.clinvarErrorMsg?.length === 0) {
                      e.stopPropagation();
                      openInNewWindow(
                        `https://www.ncbi.nlm.nih.gov/clinvar/variation/${props.vus.clinvarVariationId}`
                      );
                    }
                  }}
                />
              )}
            </div>
            <div
              className={`${styles.info} ${
                (!props.vus.clinvarClassification ||
                  props.vus.clinvarClassification?.length === 0) &&
                (!props.vus.clinvarErrorMsg ||
                  props.vus.clinvarErrorMsg?.length === 0)
                  ? styles.disabled
                  : props.vus.rsid?.length > 0 && !props.vus.rsidDbsnpVerified
                  ? styles["unverified-rsid"]
                  : ""
              }`}
            >
              {props.vus.clinvarClassification?.length > 0 ? (
                <>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Classification:</div>
                    {props.vus.clinvarClassification}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Review status:</div>
                    {props.vus.clinvarClassificationReviewStatus}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Last evaluated:</div>
                    {props.vus.clinvarClassificationLastEval}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Canonical SPDI:</div>
                    {props.vus.clinvarCanonicalSpdi}
                  </div>
                </>
              ) : props.vus.clinvarErrorMsg?.length > 0 ? (
                <div className={styles.information}>
                  <div className={styles["info-title"]}>Error message:</div>
                  {props.vus.clinvarErrorMsg}
                </div>
              ) : (
                <div className={styles.information}>
                  No Clinvar entry found based on dbSNP's RSID
                </div>
              )}
            </div>
            {props.vus.clinvarClassification?.length > 0 && (
              <div
                className={styles["clinvar-archive-link"]}
                onClick={getClinvarUpdates}
              >
                Click here to view Clinvar updates
              </div>
            )}
          </div>

          {/** DbSnp */}
          <div className={styles["info-container"]}>
            <div className={styles["external-ref-title-container"]}>
              <p className={styles["info-title"]}>dbSNP:</p>
              {props.vus.rsid?.length > 0 && (
                <Icon
                  name="external-link"
                  className={`${styles["external-link"]} ${styles.dbsnp} ${
                    props.vus.rsidDbsnpVerified
                      ? styles.active
                      : props.vus.rsid?.length > 0
                      ? styles["unverified-rsid"]
                      : styles.disabled
                  }`}
                  onClick={(e) => {
                    if (
                      props.vus.rsidDbsnpVerified ||
                      props.vus.rsid?.length > 0
                    ) {
                      e.stopPropagation();
                      openInNewWindow(
                        `https://www.ncbi.nlm.nih.gov/snp/${props.vus.rsid}`
                      );
                    }
                  }}
                />
              )}
            </div>

            <div
              className={`${styles.info} ${
                props.vus.rsidDbsnpVerified
                  ? ""
                  : props.vus.rsid?.length > 0
                  ? styles["unverified-rsid"]
                  : styles.disabled
              }`}
            >
              {props.vus.rsid?.length > 0 ? (
                <>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>
                      Is RSID verified:
                    </div>
                    {props.vus.rsidDbsnpVerified.toString()}
                  </div>
                  {props.vus.rsidDbsnpVerified ? (
                    <div className={styles.information}>
                      <div className={styles["info-title"]}>RSID:</div>
                      {props.vus.rsid}
                    </div>
                  ) : (
                    <>
                      {props.vus.rsid !== "NORSID" && (
                        <div className={styles.information}>
                          <div className={styles["info-title"]}>
                            Suggested RSID:
                          </div>
                          {props.vus.rsid}
                        </div>
                      )}
                      <div className={styles.information}>
                        <div className={styles["info-title"]}>
                          Error message:
                        </div>
                        {props.vus.rsid === "NORSID"
                          ? "No RSID found."
                          : props.vus.rsidDbsnpErrorMsgs}
                      </div>
                    </>
                  )}
                </>
              ) : (
                <div className={styles.information}>
                  No valid RSID found for this variant.
                </div>
              )}
            </div>
          </div>
        </div>

        {/** ACMG rules */}
        <div className={styles["acmg-rules"]}>
          <p className={styles["info-title"]}>ACMG rules:</p>
          <div className={styles["info-description-container"]}>
            <p className={styles["info-description"]}>
              Adding or deleting an ACMG rule requires you to fill in a
              Classification Review. You will be automatically redirected to the
              Classification Review page when trying to do so.
            </p>
            <p className={styles["info-description"]}>
              To add an ACMG rule, click on the 'plus' icon. To delete an ACMG
              rule, hover on the already-added ACMG rule and you will see a
              'bin' icon. Click on this icon.
            </p>
          </div>
          <div
            className={`${styles["acmg-rules-info"]} ${
              isAcmgEditMenuOpen ? styles["acmg-edit-open"] : ""
            }`}
          >
            <div className={styles["acmg-rules-edit"]}>
              <AcmgRulesEdit
                variantId={props.vus.id}
                variantAcmgRuleIds={props.vus.acmgRuleIds}
                allAcmgRules={props.acmgRules}
                vusService={props.vusService}
                onMenuAcmgRuleHover={(acmgRuleId?: number) =>
                  setAcmgRuleHover(acmgRuleId)
                }
                onEditIconClick={(isOpen) => setIsAcmgEditMenuOpen(isOpen)}
              />
            </div>
            {isAcmgEditMenuOpen && (
              <div className={styles["acmg-info"]}>
                <AcmgRuleInfo
                  acmgRule={props.acmgRules.find((r) => r.id === acmgRuleHover)}
                />
              </div>
            )}
          </div>
        </div>

        {/** Samples */}
        <div className={styles["samples-container"]}>
          <p className={styles["info-title"]}>Samples with this variant:</p>
          {samples.length > 0 ? (
            <>
              <div className={styles["description-container"]}>
                <p className={styles["info-description"]}>
                  Click on the sample Ids to visit the respective sample page.
                  Consequences are affected by different transcripts, i.e. HGVS,
                  according to:&nbsp;
                  <a href="https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html">
                    Ensembl
                  </a>
                  . To remove a sample from this variant, tick its checkbox and
                  click on the "Remove sample/s button".
                </p>
                <div className={styles["delete-sample"]}>
                  <Icon name="bin" fill="#008080" width={21} height={21} />
                </div>
              </div>
              <div className={styles.samples}>
                <div className={styles["samples-header"]}>
                  <span>Sample Id</span>
                  <span>HGVS</span>
                  <span>Consequence</span>
                </div>
                {samples.map((s) => (
                  <div className={styles.sample}>
                    <div className={styles["sample-info"]}>
                      <span
                        className={styles["sample-id"]}
                        onClick={() =>
                          (window.location.href = `/sample/${s.id}`)
                        }
                      >
                        <b>{s.id}</b>
                      </span>
                      <span>{s.hgvs}</span>
                      <span>
                        {s.consequence.replace("_", " ").replace("_", " ")}
                      </span>
                    </div>
                    <input
                      type="checkbox"
                      className={styles.checkbox}
                      checked={sampleIdsToRemove.some((id) => id === s.id)}
                      name={s.id.toString()}
                      onChange={() => {
                        if (sampleIdsToRemove.some((id) => id === s.id)) {
                          setSampleIdsToRemove(
                            sampleIdsToRemove.filter((id) => id !== s.id)
                          );
                        } else {
                          setSampleIdsToRemove(sampleIdsToRemove.concat(s.id));
                        }
                      }}
                    />
                  </div>
                ))}
              </div>
            </>
          ) : (
            <p className={styles["info-description"]}>
              This variant has no associated samples.
            </p>
          )}
          {sampleIdsToRemove.length > 0 && (
            <Button
              text="Remove selected sample/s"
              icon="bin"
              className={styles["remove-sample-button"]}
              onClick={() => setIsRemovingSamplesModalVisible(true)}
            />
          )}
          {notVariantSamples.length > 0 && (
            <Button
              text="Add existing samples"
              icon="add"
              className={styles["add-sample-button"]}
              onClick={() => setIsAddingSamplesModalVisible(true)}
            />
          )}
          <Button
            text="Add new sample"
            icon="add"
            className={styles["add-sample-button"]}
            onClick={() => setIsAddingNewSampleModalVisible(true)}
          />
        </div>

        {/** Phenotypes */}
        <div
          className={`${styles["phenotypes-container"]} ${
            phenotypes.length > 0 ? styles["phenotypes-populated"] : ""
          }`}
        >
          <p className={styles["info-title"]}>
            Phenotypes of the above samples:
          </p>
          {phenotypes.length > 0 ? (
            <>
              <p className={styles["info-description"]}>
                Checkout if there are any publications for this variant in
                relation to a phenotype by clicking on the book icon next to the
                phenotype.
              </p>
              <div className={styles.phenotypes}>
                {phenotypes.map((p) => (
                  <div className={styles["phenotype-container"]}>
                    <div
                      className={styles.phenotype}
                      onClick={() =>
                        openInNewWindow(
                          `https://hpo.jax.org/app/browse/term/${p.ontologyId}`
                        )
                      }
                    >
                      <b>{p.ontologyId}</b>: {p.name}
                    </div>
                    <Icon
                      name="publication"
                      className={styles["pub-icon"]}
                      onClick={() =>
                        openInNewWindow(
                          `/publication-phenotype-view/${props.vus.id}?rsid=${props.vus.rsid}&phenotype=${p.name}`
                        )
                      }
                    />
                  </div>
                ))}
              </div>
            </>
          ) : (
            <p className={styles["info-description"]}>
              The above samples have no saved phenotypes.
            </p>
          )}
        </div>
      </div>

      {showClinvarArchiveModal && (
        <Modal
          title="Clinvar updates archive"
          isClosable={true}
          onCloseIconClickCallback={() => setShowClinvarArchiveModal(false)}
        >
          {clinvarUpdates.length > 0 ? (
            <div className={styles["clinvar-updates"]}>
              <div className={styles["clinvar-updates-description"]}>
                <p>
                  Clinvar's classification was checked on every highlighted date
                  shown below.
                </p>
                <p>
                  Dates in bold indicate a change in Clinvar's germline
                  classifcation. Further details about the change can be found
                  below the calendar.
                </p>
              </div>
              <div className={styles.calendar}>
                <CalendarDisplay
                  markedDates={Array.from(
                    clinvarUpdates.map((c) => {
                      const date = c.dateChecked.split(" ")[0];

                      return {
                        date: date,
                        update:
                          datesWithUpdates?.find((d) => d === date) ?? false,
                      };
                    })
                  )}
                />
              </div>
              {clinvarUpdates
                .filter((u) => u.update !== null && u.update !== undefined)
                .map((u) => (
                  <div className={styles["clinvar-update-info-with-update"]}>
                    <p className={styles["clinvar-date-checked-container"]}>
                      <span className={styles.bullet}>{"\u25CF"}</span>
                      <span className={styles["clinvar-date-checked"]}>
                        {u.dateChecked}
                      </span>
                    </p>
                    <div className={styles["clinvar-update"]}>
                      <p>
                        <b>Last evaluated:</b> {u.update.lastEval}
                      </p>
                      <p>
                        <b>Classification:</b> {u.update.classification}
                      </p>
                      <p>
                        <b>Review status:</b> {u.update.reviewStatus}
                      </p>
                    </div>
                  </div>
                ))}
            </div>
          ) : (
            <Loader />
          )}
        </Modal>
      )}

      {isRemovingSamplesModalVisible && (
        <Modal
          title="Please Note"
          isClosable={!isRemovingSamples}
          modalContainerStyle={styles["remove-samples-modal"]}
          onCloseIconClickCallback={closeSampleRemovalModal}
        >
          <div className={styles["remove-samples-modal-content"]}>
            <p>
              Samples are deleted if they no longer have any variants. Are you
              sure you want to remove the sample/s with the following Ids:&nbsp;
              <b>{sampleIdsToRemove.join(", ")}</b> from variant&nbsp;
              <b>{props.vus.id}</b> ?
            </p>
            <div className={styles["option-btns"]}>
              <Button
                disabled={isRemovingSamples}
                text={"Yes, remove sample/s"}
                onClick={() => removeSamples()}
              />
              <Button
                disabled={isRemovingSamples}
                text="No, return to variant page"
                onClick={closeSampleRemovalModal}
              />
            </div>
          </div>
          {isRemovingSamples && <Loader />}
        </Modal>
      )}

      {isAddingSamplesModalVisible && (
        <Modal
          title={showSamplesInfoToAdd ? "Input Sample Info" : "Add Samples"}
          isClosable={!isAddingSamples}
          modalContainerStyle={`${styles["add-samples-modal"]} ${
            showSamplesInfoToAdd ? styles["show-sample-info"] : ""
          }`}
          onCloseIconClickCallback={closeAddSamplesModal}
        >
          <div className={styles["add-samples-modal-content"]}>
            {showSamplesInfoToAdd ? (
              <p>
                For each of the selected samples to be added to this variant,
                choose the genotype and input the HGVS.
              </p>
            ) : (
              <p>
                Select which of the below samples you would like to add to
                Variant&nbsp;
                <b>{props.vus.id}</b> by ticking the respective checkboxes.
              </p>
            )}

            {showSamplesInfoToAdd ? (
              <div className={styles["not-variant-samples"]}>
                {notVariantSamples
                  .filter((s) => sampleIdsToAdd.includes(s.id))
                  .map((s) => {
                    const sampleInfo =
                      sampleInfoToAdd.find(
                        (sample) => sample.sampleId === s.id
                      ) ?? null;

                    return (
                      <div className={styles["sample-to-add"]}>
                        <div className={styles.summary}>
                          <span>
                            <b>{s.id}</b>
                          </span>
                        </div>
                        <div className={styles.info}>
                          <span className={styles["info-title"]}>
                            Genotype:
                          </span>
                          <div className={styles.pills}>
                            {["Heterozygous", "Homozygous"].map((g) => {
                              return (
                                <div
                                  className={`${styles.pill} ${
                                    sampleInfo?.genotype === g
                                      ? styles["selected-pill"]
                                      : ""
                                  }`}
                                  onClick={() => {
                                    if (!isAddingSamples) {
                                      updateAddSampleGenotype(
                                        g,
                                        s.id,
                                        sampleInfo
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
                          <span className={styles["info-title"]}>HGVS:</span>
                          <Text
                            disabled={isAddingSamples}
                            onChange={(e) =>
                              updateAddSampleHgvs(
                                e.currentTarget.value,
                                s.id,
                                sampleInfo
                              )
                            }
                          />
                        </div>
                        <div className={`${styles.info} ${styles.phenotypes}`}>
                          <span className={styles["info-title"]}>
                            Phenotypes:
                          </span>
                          <SamplePhenotypeSelection
                            isDisabled={isAddingSamples}
                            isSelectingPhenotypesForVariant={true}
                            sampleService={props.sampleService}
                            onPhenotypesUpdateCallback={(phenotypes) =>
                              updateAddSamplePhenotypes(
                                phenotypes,
                                s.id,
                                sampleInfo
                              )
                            }
                          />
                        </div>
                      </div>
                    );
                  })}
              </div>
            ) : (
              <SampleTable
                sampleList={notVariantSamples.map((s) => {
                  return { sampleId: s.id, numOfVariants: s.noOfVariants };
                })}
                isClickable={false}
                showCheckboxes={true}
                onSelectedSamplesUpdate={(samplesToAdd) =>
                  setSampleIdsToAdd(samplesToAdd)
                }
              />
            )}
            <div className={styles["option-btns"]}>
              {sampleIdsToAdd.length > 0 && (
                <Button
                  disabled={
                    (showSamplesInfoToAdd &&
                      (sampleInfoToAdd.length !== sampleIdsToAdd.length ||
                        sampleInfoToAdd.filter((s) => !s.hgvs || !s.genotype)
                          .length > 0)) ||
                    isAddingSamples
                  }
                  text={
                    showSamplesInfoToAdd ? "Add samples" : "Input samples' info"
                  }
                  onClick={() => {
                    if (sampleInfoToAdd.length === 0) {
                      setShowSamplesInfoToAdd(true);
                    } else {
                      addSamples();
                    }
                  }}
                />
              )}
              <Button
                text="Return to variant page"
                onClick={closeAddSamplesModal}
                disabled={isAddingSamples}
              />
            </div>
          </div>
          {isAddingSamples && <Loader />}
        </Modal>
      )}

      {isAddingNewSampleModalVisible && (
        <Modal
          title={"Input Sample Info"}
          isClosable={!isAddingSamples}
          modalContainerStyle={`${styles["add-new-sample-modal"]} ${styles["show-sample-info"]}`}
          onCloseIconClickCallback={closeAddNewSampleModal}
        >
          <div className={styles["add-new-sample-modal-content"]}>
            <p>Input the sample Id, select its genotype and input its HGVS.</p>

            <div className={styles["not-variant-samples"]}>
              <div className={styles["sample-to-add"]}>
                <div>
                  <div className={styles.info}>
                    <span className={styles["info-title"]}>Id:</span>
                    <Text
                      disabled={isAddingSamples}
                      onChange={(e) => {
                        setSampleInfoToAdd([
                          {
                            ...sampleInfoToAdd[0],
                            sampleId: e.currentTarget.value,
                          },
                        ]);
                        setShowNewSampleIdError(
                          e.currentTarget.value.trim().length > 0 &&
                            (samples
                              .map((s) => s.id)
                              .includes(e.currentTarget.value) ||
                              notVariantSamples
                                .map((s) => s.id)
                                .includes(e.currentTarget.value))
                        );
                      }}
                    />
                  </div>
                  {showNewSampleIdError && (
                    <p className={styles.error}>This sample already exists!</p>
                  )}
                </div>
                <div className={styles.info}>
                  <span className={styles["info-title"]}>Genotype:</span>
                  <div className={styles.pills}>
                    {["Heterozygous", "Homozygous"].map((g) => {
                      return (
                        <div
                          className={`${styles.pill} ${
                            sampleInfoToAdd[0]?.genotype === g
                              ? styles["selected-pill"]
                              : ""
                          }`}
                          onClick={() => {
                            if (!isAddingSamples) {
                              setSampleInfoToAdd([
                                {
                                  ...sampleInfoToAdd[0],
                                  genotype: g,
                                },
                              ]);
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
                  <span className={styles["info-title"]}>HGVS:</span>
                  <Text
                    disabled={isAddingSamples}
                    onChange={(e) =>
                      setSampleInfoToAdd([
                        {
                          ...sampleInfoToAdd[0],
                          hgvs: e.currentTarget.value,
                        },
                      ])
                    }
                  />
                </div>
                <div className={`${styles.info} ${styles.phenotypes}`}>
                  <span className={styles["info-title"]}>Phenotypes:</span>
                  <SamplePhenotypeSelection
                    isDisabled={isAddingSamples}
                    isSelectingPhenotypesForVariant={true}
                    sampleService={props.sampleService}
                    onPhenotypesUpdateCallback={(phenotypes) =>
                      setSampleInfoToAdd([
                        {
                          ...sampleInfoToAdd[0],
                          phenotypes: phenotypes,
                        },
                      ])
                    }
                  />
                </div>
              </div>
            </div>

            <div className={styles["option-btns"]}>
              <Button
                disabled={
                  sampleInfoToAdd.length !== 1 ||
                  sampleInfoToAdd.filter(
                    (s) => !s.sampleId || !s.hgvs || !s.genotype
                  ).length > 0 ||
                  showNewSampleIdError ||
                  isAddingSamples
                }
                text={"Add new sample"}
                onClick={() => {
                  addNewSample();
                }}
              />
              <Button
                text="Return to variant page"
                onClick={closeAddNewSampleModal}
                disabled={isAddingSamples}
              />
            </div>
          </div>
          {isAddingSamples && <Loader />}
        </Modal>
      )}
    </div>
  );

  function getClinvarUpdates() {
    setShowClinvarArchiveModal(true);

    if (clinvarUpdates.length === 0) {
      props.vusService
        .getClinvarUpdates({ clinvarId: props.vus.clinvarId })
        .then((res) => {
          if (res.isSuccess) {
            setClinvarUpdates(res.clinvarUpdates);
            if (res.datesWithUpdates) setDatesWithUpdates(res.datesWithUpdates);
          }
        });
    }
  }

  function updateAddSampleGenotype(
    genotype: string,
    sampleId: string,
    sampleInfo?: ISampleToAddInfo
  ) {
    if (sampleInfo) {
      let samplesInfoToAddUpdated = [];

      sampleInfoToAdd.forEach((info) => {
        if (info.sampleId === sampleId) {
          samplesInfoToAddUpdated = samplesInfoToAddUpdated.concat({
            ...info,
            genotype: genotype,
          });
        } else {
          samplesInfoToAddUpdated = samplesInfoToAddUpdated.concat(info);
        }
      });

      setSampleInfoToAdd(samplesInfoToAddUpdated);
    } else {
      setSampleInfoToAdd(
        sampleInfoToAdd.concat({
          sampleId: sampleId,
          genotype: genotype,
        })
      );
    }
  }

  function updateAddSampleHgvs(
    hgvs: string,
    sampleId: string,
    sampleInfo?: ISampleToAddInfo
  ) {
    if (sampleInfo) {
      let samplesInfoToAddUpdated = [];

      sampleInfoToAdd.forEach((info) => {
        if (info.sampleId === sampleId) {
          samplesInfoToAddUpdated = samplesInfoToAddUpdated.concat({
            ...info,
            hgvs: hgvs,
          });
        } else {
          samplesInfoToAddUpdated = samplesInfoToAddUpdated.concat(info);
        }
      });

      setSampleInfoToAdd(samplesInfoToAddUpdated);
    } else {
      setSampleInfoToAdd(
        sampleInfoToAdd.concat({
          sampleId: sampleId,
          hgvs: hgvs,
        })
      );
    }
  }

  function updateAddSamplePhenotypes(
    phenotypes: IPhenotype[],
    sampleId: string,
    sampleInfo?: ISampleToAddInfo
  ) {
    if (sampleInfo) {
      let samplesInfoToAddUpdated = [];

      sampleInfoToAdd.forEach((info) => {
        if (info.sampleId === sampleId) {
          samplesInfoToAddUpdated = samplesInfoToAddUpdated.concat({
            ...info,
            phenotypes: phenotypes,
          });
        } else {
          samplesInfoToAddUpdated = samplesInfoToAddUpdated.concat(info);
        }
      });

      setSampleInfoToAdd(samplesInfoToAddUpdated);
    } else {
      setSampleInfoToAdd(
        sampleInfoToAdd.concat({
          sampleId: sampleId,
          phenotypes: phenotypes,
        })
      );
    }
  }

  function closeAddSamplesModal() {
    setIsAddingSamplesModalVisible(false);
    setSampleIdsToAdd([]);
    setShowSamplesInfoToAdd(false);
    setSampleInfoToAdd([]);
  }

  function closeAddNewSampleModal() {
    setIsAddingNewSampleModalVisible(false);
    setSampleInfoToAdd([]);
  }

  function addSamples() {
    setIsAddingSamples(true);

    props.vusService
      .addSamples({
        variantId: props.vus.id,
        samplesToAdd: sampleInfoToAdd,
      })
      .then((res) => {
        if (res.isSuccess) {
          setSamples(res.updatedSamples);
          setNotVariantSamples(res.updatedNotVariantSamples);
          setPhenotypes(res.updatedPhenotypes);
          closeAddSamplesModal();
        }
        setIsAddingSamples(false);
      });
  }

  function addNewSample() {
    setIsAddingSamples(true);

    props.vusService
      .addNewSample({
        variantId: props.vus.id,
        sampleToAdd: sampleInfoToAdd[0],
      })
      .then((res) => {
        if (res.isSuccess) {
          setSamples(res.updatedSamples);
          setNotVariantSamples(res.updatedNotVariantSamples);
          setPhenotypes(res.updatedPhenotypes);
          closeAddNewSampleModal();
        }
        setIsAddingSamples(false);
      });
  }

  function closeSampleRemovalModal() {
    setIsRemovingSamples(false);
    setIsRemovingSamplesModalVisible(false);
    setSampleIdsToRemove([]);
  }

  function removeSamples() {
    setIsRemovingSamples(true);

    props.vusService
      .removeSamples({
        variantId: props.vus.id,
        sampleIdsToRemove: sampleIdsToRemove,
      })
      .then((res) => {
        if (res.isSuccess) {
          setSamples(res.updatedSamples);
          setNotVariantSamples(res.updatedNotVariantSamples);
          setPhenotypes(res.updatedPhenotypes);
          closeSampleRemovalModal();
        }
        setIsRemovingSamples(false);
      });
  }
};

export default VusInfo;
