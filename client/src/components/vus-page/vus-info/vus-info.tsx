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
import DebouncedInput from "../../shared/table/debounced-input/debounced-input";

type VusInfoProps = {
  vus: IVus;
  acmgRules: IAcmgRule[];
  vusService?: VusService;
  sampleService?: SampleService;
  onRsidUpdateCallback?: (numOfPublications: number) => void;
};

const VusInfo: React.FunctionComponent<VusInfoProps> = (
  props: VusInfoProps
) => {
  const [filterValue, setFilterValue] = useState("");

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

  const [updatedRsid, setUpdatedRsid] = useState(undefined);
  const [updatedRsidErrorMsg, setUpdatedRsidErrorMsg] = useState("");
  const [showUpdateRsidModal, setShowUpdateRsidModal] = useState(false);
  const [showUpdateRsidConfirmationModal, setShowUpdateRsidConfirmationModal] =
    useState(false);
  const [isUpdatingRsid, setIsUpdatingRsid] = useState(false);

  const [rsidInfo, setRsidInfo] = useState({
    rsid: props.vus.rsid,
    rsidDbsnpVerified: props.vus.rsidDbsnpVerified,
    rsidDbsnpErrorMsgs: props.vus.rsidDbsnpErrorMsgs,
  });
  const [clinvarInfo, setClinvarInfo] = useState({
    clinvarId: props.vus.clinvarId,
    clinvarVariationId: props.vus.clinvarVariationId,
    clinvarCanonicalSpdi: props.vus.clinvarCanonicalSpdi,
    clinvarClassification: props.vus.clinvarClassification,
    clinvarClassificationReviewStatus:
      props.vus.clinvarClassificationReviewStatus,
    clinvarClassificationLastEval: props.vus.clinvarClassificationLastEval,
    clinvarErrorMsg: props.vus.clinvarErrorMsg,
  });

  const rsidEditIcon = (
    <Icon
      width={16}
      height={16}
      name="edit"
      className={styles["edit-rsid-icon"]}
      onClick={() => setShowUpdateRsidModal(true)}
    />
  );

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

        {/** External References */}
        <div className={styles["external-ref"]}>
          {/** Clinvar */}
          <div className={styles["info-container"]}>
            <div className={styles["external-ref-title-container"]}>
              <p className={styles["info-title"]}>
                {`ClinVar${
                  clinvarInfo.clinvarClassification?.length > 0 &&
                  !rsidInfo.rsidDbsnpVerified
                    ? " of suggested dbSNP RSID"
                    : ""
                }:`}
              </p>
              {(clinvarInfo.clinvarClassification?.length > 0 ||
                clinvarInfo.clinvarErrorMsg?.length > 0) && (
                <Icon
                  width={18}
                  height={18}
                  name="external-link"
                  className={`${styles["external-link"]} ${styles.clinvar} ${
                    clinvarInfo.clinvarClassification?.length === 0
                      ? styles.disabled
                      : (rsidInfo.rsid?.length > 0 &&
                          !rsidInfo.rsidDbsnpVerified) ||
                        clinvarInfo.clinvarErrorMsg.length > 0
                      ? styles["unverified-rsid"]
                      : styles.active
                  }`}
                  onClick={(e) => {
                    if (clinvarInfo.clinvarVariationId.length > 0) {
                      e.stopPropagation();
                      openInNewWindow(
                        `https://www.ncbi.nlm.nih.gov/clinvar/variation/${clinvarInfo.clinvarVariationId}`
                      );
                    }
                  }}
                />
              )}
            </div>
            <div
              className={`${styles.info} ${
                (!clinvarInfo.clinvarClassification ||
                  clinvarInfo.clinvarClassification?.length === 0) &&
                (!clinvarInfo.clinvarErrorMsg ||
                  clinvarInfo.clinvarErrorMsg?.length === 0)
                  ? styles.disabled
                  : (rsidInfo.rsid?.length > 0 &&
                      !rsidInfo.rsidDbsnpVerified) ||
                    clinvarInfo.clinvarErrorMsg.length > 0
                  ? styles["unverified-rsid"]
                  : ""
              }`}
            >
              {clinvarInfo.clinvarClassification?.length > 0 && (
                <>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Classification:</div>
                    {clinvarInfo.clinvarClassification}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Review status:</div>
                    {clinvarInfo.clinvarClassificationReviewStatus}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Last evaluated:</div>
                    {clinvarInfo.clinvarClassificationLastEval}
                  </div>
                  <div className={styles.information}>
                    <div className={styles["info-title"]}>Canonical SPDI:</div>
                    <span style={{ overflowX: "auto" }}>
                      {clinvarInfo.clinvarCanonicalSpdi}
                    </span>
                  </div>
                </>
              )}

              {clinvarInfo.clinvarErrorMsg?.length > 0 && (
                <div className={styles.information}>
                  <div className={styles["info-title"]}>Error message:</div>
                  {clinvarInfo.clinvarErrorMsg}
                </div>
              )}

              {(!clinvarInfo?.clinvarClassification ||
                clinvarInfo?.clinvarClassification?.length === 0) && (
                <div className={styles.information}>
                  No Clinvar entry found based on dbSNP's RSID
                </div>
              )}
            </div>
            {clinvarInfo.clinvarClassification?.length > 0 && (
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
              {rsidInfo.rsid?.length > 0 && (
                <Icon
                  width={18}
                  height={18}
                  name="external-link"
                  className={`${styles["external-link"]} ${styles.dbsnp} ${
                    rsidInfo.rsidDbsnpVerified
                      ? styles.active
                      : rsidInfo.rsid?.length > 0
                      ? styles["unverified-rsid"]
                      : styles.disabled
                  }`}
                  onClick={(e) => {
                    if (
                      rsidInfo.rsidDbsnpVerified ||
                      rsidInfo.rsid?.length > 0
                    ) {
                      e.stopPropagation();
                      openInNewWindow(
                        `https://www.ncbi.nlm.nih.gov/snp/${rsidInfo.rsid}`
                      );
                    }
                  }}
                />
              )}
            </div>

            <div
              className={`${styles.info} ${
                rsidInfo.rsidDbsnpVerified
                  ? ""
                  : rsidInfo.rsid?.length > 0
                  ? styles["unverified-rsid"]
                  : styles.disabled
              }`}
            >
              {rsidInfo.rsid?.length > 0 ? (
                <>
                  {rsidInfo.rsidDbsnpVerified ? (
                    <div className={styles["rsid-information"]}>
                      <div className={styles.information}>
                        <div className={styles["info-title"]}>RSID:</div>
                        {rsidInfo.rsid}
                      </div>
                      {rsidEditIcon}
                    </div>
                  ) : (
                    <>
                      {rsidInfo.rsid !== "NORSID" && (
                        <div className={styles.information}>
                          <div className={styles["info-title"]}>
                            Suggested RSID:
                          </div>
                          {rsidInfo.rsid}
                        </div>
                      )}
                      <div className={styles.information}>
                        <div className={styles["info-title"]}>
                          Error message/s:
                        </div>
                        {rsidInfo.rsid === "NORSID" ? (
                          "No RSID found."
                        ) : (
                          <div className={styles["rsid-errors"]}>
                            {rsidInfo.rsidDbsnpErrorMsgs
                              .split("|| ")
                              .map((msg) => (
                                <div className={styles["rsid-error"]}>
                                  <div className={styles.bullet}>
                                    {"\u25CF"}
                                  </div>
                                  <span>{msg}</span>
                                </div>
                              ))}
                          </div>
                        )}
                      </div>
                    </>
                  )}
                </>
              ) : (
                <div className={styles["rsid-information"]}>
                  <p>No valid RSID found for this variant.</p>
                  {rsidEditIcon}
                </div>
              )}
            </div>

            {rsidInfo.rsidDbsnpErrorMsgs?.length > 0 && (
              <>
                <Button
                  text="Save Suggested RSID"
                  className={styles["rsid-update-btn"]}
                  onClick={() => {
                    setUpdatedRsid(props.vus.rsid);
                    setShowUpdateRsidConfirmationModal(true);
                  }}
                />
                <Button
                  text="Update RSID"
                  className={styles["rsid-update-btn"]}
                  onClick={() => setShowUpdateRsidModal(true)}
                />
              </>
            )}
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
          <p className={styles["info-title"]}>Variant samples:</p>
          {samples.length > 0 ? (
            <>
              <p className={styles["info-description"]}>
                Click on the sample Ids to visit the respective sample page.
                Consequences are affected by different transcripts, i.e. HGVS,
                according to:&nbsp;
                <a href="https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html">
                  Ensembl
                </a>
                . To remove a sample from this variant, tick its checkbox and
                click on the "Remove selected sample/s" button at the bottom of
                the sample list. You can search through the samples using sample
                ids, HGVS or consequences.
              </p>
              <DebouncedInput
                onChange={(val) => setFilterValue(val.toString())}
                placeholder={`Search samples ...`}
                type="text"
                value={filterValue}
                className={styles.input}
              />
              <div className={styles.samples}>
                <div className={styles["samples-header"]}>
                  <span>Sample Id</span>
                  <span>HGVS</span>
                  <span>Consequence</span>
                  <div className={styles["delete-sample"]}>
                    <Icon name="bin" fill="#008080" width={21} height={21} />
                  </div>
                </div>
                {(filterValue.length > 0
                  ? samples.filter(
                      (s) =>
                        containsFilterValue(s.id) ||
                        containsFilterValue(s.hgvs) ||
                        containsFilterValue(s.consequence)
                    )
                  : samples
                ).map((s) => (
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
                        {s.consequence?.replace("_", " ").replace("_", " ")}
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
                    {/** Allow user to look up publications related to the phenotype, only if the vus has at least 1 HGVS or has a valid RSID  */}
                    {((rsidInfo.rsid &&
                      rsidInfo.rsidDbsnpErrorMsgs?.length === 0) ||
                      samples.find((s) => s.hgvs?.length > 0)) && (
                      <Icon
                        name="publication"
                        className={styles["pub-icon"]}
                        onClick={() =>
                          openInNewWindow(
                            `/publication-phenotype-view/${props.vus.id}?rsid=${rsidInfo.rsid}&phenotype=${p.name}`
                          )
                        }
                      />
                    )}
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
                        sampleInfoToAdd.filter((s) => !s.genotype).length >
                          0)) ||
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
                  sampleInfoToAdd.filter((s) => !s.sampleId || !s.genotype)
                    .length > 0 ||
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

      {showUpdateRsidModal && (
        <Modal
          modalContainerStyle={styles["update-rsid-confirmation-modal"]}
          title="Update RSID"
        >
          <div className={styles["update-rsid-confirmation-modal-content"]}>
            <div className={styles["update-rsid"]}>
              <p>RSID:</p>
              <div className={styles["rsid-input"]}>
                <Text
                  value={updatedRsid}
                  autoFocus={true}
                  maxLength={15}
                  onChange={(e) => {
                    setUpdatedRsid(e.currentTarget.value);
                    setUpdatedRsidErrorMsg(
                      (e.currentTarget.value.length === 1 &&
                        !e.currentTarget.value.startsWith("r")) ||
                        (e.currentTarget.value.length > 1 &&
                          !e.currentTarget.value.match(/^rs\d{1,13}$/))
                        ? "RSID needs to start with 'rs' followed by a maximum of 13 digits"
                        : ""
                    );
                  }}
                />
                {updatedRsidErrorMsg.length > 0 && (
                  <p className={styles.error}>{updatedRsidErrorMsg}</p>
                )}
              </div>
            </div>
            <div className={styles["option-btns"]}>
              <Button
                text="Update RSID"
                onClick={() => {
                  setShowUpdateRsidModal(false);
                  setShowUpdateRsidConfirmationModal(true);
                }}
                className={styles["delete-btn"]}
                disabled={
                  isUpdatingRsid ||
                  !updatedRsid ||
                  updatedRsid.trim().length === 0 ||
                  updatedRsidErrorMsg.length > 0
                }
              />
              <Button
                text="Cancel update"
                onClick={() => {
                  setIsUpdatingRsid(false);
                  setShowUpdateRsidModal(false);
                  setUpdatedRsid(undefined);
                  setUpdatedRsidErrorMsg("");
                }}
                disabled={isUpdatingRsid}
              />
            </div>
          </div>
          {isUpdatingRsid && <Loader />}
        </Modal>
      )}

      {showUpdateRsidConfirmationModal && (
        <Modal modalContainerStyle={styles["update-rsid-confirmation-modal"]}>
          <div className={styles["update-rsid-confirmation-modal-content"]}>
            <p>
              Are you sure you want to update Variant&nbsp;
              <b>{props.vus.id}</b>'s RSID to <b>{updatedRsid}</b> ?
            </p>
            <div className={styles["option-btns"]}>
              <Button
                text="Yes, update it"
                onClick={saveSuggestedRsid}
                className={styles["delete-btn"]}
                disabled={isUpdatingRsid}
              />
              <Button
                text="No, return to variant page"
                onClick={() => {
                  setIsUpdatingRsid(false);
                  setShowUpdateRsidConfirmationModal(false);
                  setUpdatedRsid(undefined);
                  setUpdatedRsidErrorMsg("");
                }}
                disabled={isUpdatingRsid}
              />
            </div>
          </div>
          {isUpdatingRsid && <Loader />}
        </Modal>
      )}
    </div>
  );

  function getClinvarUpdates() {
    setShowClinvarArchiveModal(true);

    if (clinvarUpdates.length === 0) {
      props.vusService
        .getClinvarUpdates({ clinvarId: clinvarInfo.clinvarId })
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

  function saveSuggestedRsid() {
    setIsUpdatingRsid(true);
    props.vusService
      .updateRsid({
        variantId: props.vus.id,
        newRsid: updatedRsid,
      })
      .then((res) => {
        if (res.isSuccess) {
          setRsidInfo({
            rsid: res.updatedExternalRefData.rsid,
            rsidDbsnpErrorMsgs: res.updatedExternalRefData.rsidDbsnpErrorMsgs,
            rsidDbsnpVerified: res.updatedExternalRefData.rsidDbsnpVerified,
          });

          setClinvarInfo({
            clinvarId: res.updatedExternalRefData.clinvarId,
            clinvarVariationId: res.updatedExternalRefData.clinvarVariationId,
            clinvarCanonicalSpdi:
              res.updatedExternalRefData.clinvarCanonicalSpdi,
            clinvarClassification:
              res.updatedExternalRefData.clinvarClassification,
            clinvarClassificationLastEval:
              res.updatedExternalRefData.clinvarClassificationLastEval,
            clinvarErrorMsg: res.updatedExternalRefData.clinvarErrorMsg,
            clinvarClassificationReviewStatus:
              res.updatedExternalRefData.clinvarClassificationReviewStatus,
          });

          setClinvarUpdates([]);

          if (res.updatedExternalRefData.clinvarId) {
            props.vusService
              .getClinvarUpdates({
                clinvarId: res.updatedExternalRefData.clinvarId,
              })
              .then((res) => {
                if (res.isSuccess) {
                  setClinvarUpdates(res.clinvarUpdates);
                  if (res.datesWithUpdates)
                    setDatesWithUpdates(res.datesWithUpdates);
                }
              });
          }

          props.onRsidUpdateCallback &&
            props.onRsidUpdateCallback(
              res.updatedExternalRefData.numOfPublications
            );
        }
        setShowUpdateRsidConfirmationModal(false);
        setIsUpdatingRsid(false);
        setUpdatedRsid(undefined);
        setUpdatedRsidErrorMsg("");
      });
  }

  function containsFilterValue(val: string): boolean {
    if (val && val.length > 0) {
      return val?.toLowerCase().includes(filterValue?.toLowerCase());
    }
    return false;
  }
};

export default VusInfo;
