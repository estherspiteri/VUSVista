import React, { useState } from "react";
import styles from "./vus-upload-page.module.scss";
import Text from "./../../atoms/text/text";
import VusUploadField from "./vus-upload-field/vus-upload-field";
import Button from "../../atoms/button/button";
import SamplePhenotypeSelection from "../sample-page/sample-phenotype-selection/sample-phenotype-selection";
import { IHPOTerm } from "../../services/sample/sample.dto";
import { SampleService } from "../../services/sample/sample.service";
import Icon from "../../atoms/icons/icon";
import { VusService } from "../../services/vus/vus.service";
import Modal from "../../atoms/modal/modal";
import { IAcmgRuleUpload, IVusUpload } from "../../models/vus-upload.model";
import AcmgRuleInfo from "../sample-page/acmg-rule-info/acmg-rule-info";
import AcmgRulesEdit from "../sample-page/acmg-rules-edit/acmg-rules-edit";
import { IAcmgRule } from "../../models/acmg-rule.model";
import Loader from "../../atoms/loader/loader";

type VusUploadPageProps = {
  acmgRules: IAcmgRule[];
  sampleService?: SampleService;
  vusService?: VusService;
};

const VusUploadPage: React.FunctionComponent<VusUploadPageProps> = (
  props: VusUploadPageProps
) => {
  const [acmgRuleHover, setAcmgRuleHover] = useState<number | undefined>(
    undefined
  );

  const [isLoading, setIsLoading] = useState(false);

  const [chromosome, setChromosome] = useState<string | undefined>(undefined);
  const [chromosomeErrorMsg, setChromosomeErrorMsg] = useState("");

  const [chromosomePosition, setChromosomePosition] = useState<
    number | undefined
  >(undefined);
  const [chromosomePositionErrorMsg, setChromosomePositionErrorMsg] =
    useState("");

  const [refAllele, setRefAllele] = useState<string | undefined>("");
  const [refAlleleErrorMsg, setRefAlleleErrorMsg] = useState("");

  const [altAllele, setAltAllele] = useState<string | undefined>("");
  const [altAlleleErrorMsg, setAltAlleleErrorMsg] = useState("");

  const [genotype, setGenotype] = useState<string | undefined>(undefined);
  const [genotypeErrorMsg, setGenotypeErrorMsg] = useState("");

  const [geneInput, setGeneInput] = useState<string | undefined>(undefined);
  const [isValidatingGeneInput, setIsValidatingGeneInput] = useState(false);
  const [geneInputErrorMsg, setGeneInputErrorMsg] = useState("");
  const [geneId, setGeneId] = useState<number | undefined>(undefined);

  const [classification, setClassification] = useState<string | undefined>(
    undefined
  );
  const [classificationErrorMsg, setClassificationErrorMsg] = useState("");

  const [type, setType] = useState<string | undefined>(undefined);
  const [typeErrorMsg, setTypeErrorMsg] = useState("");

  const [samples, setSamples] = useState<string[]>([]);
  const [samplesErrorMsg, setSamplesErrorMsg] = useState("");

  const [phenotypes, setPhenotypes] = useState<IHPOTerm[]>(undefined);

  const [acmgRules, setAcmgRules] = useState<IAcmgRuleUpload[]>(undefined);

  const [isSummaryModalOpen, setIsSummaryModelOpen] = useState(false);

  //validity of the sections
  const isChrValid =
    chromosome !== undefined &&
    chromosomeErrorMsg.length === 0 &&
    chromosomePosition !== undefined &&
    chromosomePositionErrorMsg.length === 0;

  const areAllelesValid =
    refAllele.length > 0 &&
    refAlleleErrorMsg.length === 0 &&
    altAllele.length > 0 &&
    altAlleleErrorMsg.length === 0;

  const isGenotypeValid =
    genotype !== undefined && genotypeErrorMsg.length === 0;

  const isGeneValid = geneId !== undefined;

  const isClassificationValid =
    classification !== undefined && classificationErrorMsg.length === 0;

  const isTypeValid = type !== undefined && typeErrorMsg.length === 0;

  const areSamplesValid = samples.length > 0;

  if (isLoading) {
    return <Loader />;
  } else {
    return (
      <div className={styles["vus-upload-page-container"]}>
        <div className={styles.title}>VUS Upload</div>
        <div className={styles.description}>
          <p>
            Click on each section to fill in the variant details. A correctly
            filled-in section is marked with a tick on the right-hand side.{" "}
          </p>
        </div>

        <div className={styles["fields-wrapper"]}>
          {/** Chromsome & Chromsome Position*/}
          <VusUploadField
            title="Chromosome"
            isOpenByDefault={true}
            showCheckMark={isChrValid}
          >
            <div className={styles["field-content-container"]}>
              <div className={styles["field-content"]}>
                <span>Select chromosome: </span>
                <div className={styles.chromosomes}>
                  {[
                    "1",
                    "2",
                    "3",
                    "4",
                    "5",
                    "6",
                    "7",
                    "8",
                    "9",
                    "10",
                    "11",
                    "12",
                    "13",
                    "14",
                    "15",
                    "16",
                    "17",
                    "18",
                    "19",
                    "20",
                    "21",
                    "22",
                    "X",
                    "Y",
                  ].map((c) => (
                    <div
                      className={`${styles.chromosome} ${
                        c === chromosome ? styles["selected-chromosome"] : ""
                      }`}
                      onClick={() => setChromosome(c)}
                    >
                      {c}
                    </div>
                  ))}
                </div>
              </div>

              <div className={styles["field-content"]}>
                <span>Chromosome Position: </span>
                <Text
                  value={chromosomePosition}
                  placeholder="0"
                  type="number"
                  min={0}
                  className={styles["chromosome-position"]}
                  onChange={(e) =>
                    setChromosomePosition(parseInt(e.currentTarget.value))
                  }
                  validationCallback={checkChrPositionValidity}
                  errorMsg={chromosomePositionErrorMsg}
                />
              </div>
            </div>
          </VusUploadField>

          {/** Reference & Alt Alleles */}
          {/**TODO: check on attempt to submit if allele consists of just GACT*/}
          <VusUploadField title="Alleles" showCheckMark={areAllelesValid}>
            <div className={styles["allele-wrapper"]}>
              <div className={styles["field-content"]}>
                <span>Reference allele:</span>
                <Text
                  value={refAllele}
                  onChange={(e) => setRefAllele(e.currentTarget.value)}
                  validationCallback={(val) => checkAlleleValidity(val, true)}
                  errorMsg={refAlleleErrorMsg}
                />
              </div>
              <div className={styles["field-content"]}>
                <span className={styles["alt-allele-description"]}>
                  Alternate allele ('Input '/' if there is no alternate'):
                </span>
                <Text
                  value={altAllele}
                  onChange={(e) => setAltAllele(e.currentTarget.value)}
                  validationCallback={(val) => checkAlleleValidity(val, false)}
                  errorMsg={altAlleleErrorMsg}
                />
              </div>
            </div>
          </VusUploadField>

          {/** Zygosity (Genotype) */}
          <VusUploadField title="Zygosity" showCheckMark={isGenotypeValid}>
            <div className={styles["field-content"]}>
              <span>Select the zygosity:</span>
              <div className={styles.pills}>
                {["Heterozygous", "Homozygous"].map((g) => {
                  return (
                    <div
                      className={`${styles.pill} ${
                        genotype === g ? styles["selected-pill"] : ""
                      }`}
                      onClick={() => setGenotype(g)}
                    >
                      {g}
                    </div>
                  );
                })}
              </div>
            </div>
          </VusUploadField>

          {/** Gene */}
          <VusUploadField title="Gene" showCheckMark={isGeneValid}>
            <div className={styles["field-content"]}>
              {geneId === undefined ? (
                <>
                  <span>Type the gene and wait for it to be validated:</span>
                  <Text
                    value={geneInput}
                    disabled={isValidatingGeneInput}
                    onChange={(e) => setGeneInput(e.currentTarget.value)}
                    validationCallback={checkGeneValidity}
                    errorMsg={geneInputErrorMsg}
                  />
                </>
              ) : (
                <>
                  <span>The variant is found in the below gene:</span>
                  <div className={`${styles["selected-gene"]} ${styles.pill}`}>
                    <span>{geneInput.toUpperCase()}</span>
                    <Icon
                      name="close"
                      onClick={() => removeGene()}
                      stroke="#008080"
                    />
                  </div>
                </>
              )}
            </div>
          </VusUploadField>

          {/** Classification */}
          <VusUploadField
            title="Classification"
            showCheckMark={isClassificationValid}
          >
            <div className={styles["field-content"]}>
              <span>Select the classification:</span>
              <div className={styles.pills}>
                {["VUS", "Uncertain significance", "Unclassified"].map((c) => {
                  return (
                    <div
                      className={`${styles.pill} ${
                        classification === c ? styles["selected-pill"] : ""
                      }`}
                      onClick={() => setClassification(c)}
                    >
                      {c}
                    </div>
                  );
                })}
              </div>
            </div>
          </VusUploadField>

          {/** Type */}
          <VusUploadField title="Type" showCheckMark={isTypeValid}>
            <div className={styles["field-content"]}>
              <span>Select the type:</span>
              <div className={styles.pills}>
                {["SNV", "MNV", "INDEL"].map((c) => {
                  return (
                    <div
                      className={`${styles.pill} ${
                        type === c ? styles["selected-pill"] : ""
                      }`}
                      onClick={() => setType(c)}
                    >
                      {c}
                    </div>
                  );
                })}
              </div>
            </div>
          </VusUploadField>

          {/** Samples */}
          <VusUploadField title="Samples" showCheckMark={areSamplesValid}>
            <div className={styles["field-content-container"]}>
              {/** Sample Ids */}
              <div className={styles["field-content"]}>
                <span>
                  List the samples that have this variant below. Press 'Enter'
                  to add the sample to the list.
                </span>
                <div className={styles.samples}>
                  <input
                    id="sample-input"
                    className={styles["sample-input"]}
                    type="text"
                    placeholder="Type in a sample id . . ."
                    onKeyDown={(e) => {
                      if (e.key === "Enter") addSample(e.currentTarget.value);
                    }}
                  />
                  <div className={styles["samples-selected"]}>
                    {samplesErrorMsg !== undefined &&
                      samples.length > 0 &&
                      samples.map((s) => {
                        return (
                          <p className={styles.sample}>
                            <Icon
                              name="close"
                              onClick={() => removeSample(s)}
                              stroke="#008080"
                            />
                            <span className={styles["selected-sample"]}>
                              {s}
                            </span>
                          </p>
                        );
                      })}
                  </div>
                </div>
              </div>

              {/** Phenotypes */}
              <div className={styles["field-content"]}>
                <span>Input the samples' phenotypes</span>
                <SamplePhenotypeSelection
                  isSelectingPhenotypesForVariant={true}
                  sampleService={props.sampleService}
                  onPhenotypesUpdateCallback={(phenotypes) =>
                    setPhenotypes(phenotypes)
                  }
                />
              </div>
            </div>
          </VusUploadField>

          {/** ACMG Rules */}
          <VusUploadField title="ACMG Rules">
            {/** ACMG Rules */}
            <div className={styles["field-content"]}>
              <span>Choose the ACMG rules that this variant has</span>
              <div className={styles["acmg-rules-edit"]}>
                <AcmgRulesEdit
                  allAcmgRules={props.acmgRules}
                  isAcmgMenuClosable={false}
                  onMenuAcmgRuleHover={(acmgRuleId?: number) =>
                    setAcmgRuleHover(acmgRuleId)
                  }
                  onAcmgRulesSelectionUpdate={(rules) => {
                    setAcmgRules(
                      props.acmgRules
                        .filter((r) => rules.includes(r.id))
                        .map((r) => {
                          return { id: r.id, name: r.name };
                        })
                    );
                  }}
                />
              </div>
              <div className={styles["acmg-rules-info"]}>
                <AcmgRuleInfo
                  acmgRule={props.acmgRules?.find(
                    (r) => r.id === acmgRuleHover
                  )}
                />
              </div>
            </div>
          </VusUploadField>
        </div>

        <Button
          className={styles.btn}
          text="Save variant"
          disabled={
            !isChrValid ||
            !areAllelesValid ||
            !isGenotypeValid ||
            !isGeneValid ||
            !isTypeValid ||
            !isClassificationValid ||
            !areSamplesValid
          }
          onClick={() => setIsSummaryModelOpen(true)}
        />

        {/** Summary Modal */}
        {isSummaryModalOpen && (
          <Modal
            title="Variant Summary"
            modalContainerStyle={styles["summary-modal"]}
          >
            <div className={styles["summary-modal-content"]}>
              <div className={styles.summary}>
                <div className={styles["summary-field"]}>
                  <p className={styles["selection-name"]}>Locus</p>
                  <span>
                    <b>
                      chr{chromosome}:{chromosomePosition}
                    </b>
                  </span>
                </div>
                <div className={styles["summary-field"]}>
                  <p className={styles["selection-name"]}>Reference allele</p>
                  <span>
                    <b>{refAllele.toUpperCase()}</b>
                  </span>
                </div>
                <div className={styles["summary-field"]}>
                  <p className={styles["selection-name"]}>Alternate allele</p>
                  <span>
                    <b>{altAllele.toUpperCase()}</b>
                  </span>
                </div>
                <div className={styles["summary-field"]}>
                  <p className={styles["selection-name"]}>Gene</p>
                  <span>
                    <b>{geneInput}</b>
                  </span>
                </div>
                <div className={styles["summary-field"]}>
                  <p className={styles["selection-name"]}>Classification</p>
                  <span>
                    <b>{classification}</b>
                  </span>
                </div>

                <div className={styles["summary-field"]}>
                  <p className={styles["selection-name"]}>Sample Ids</p>
                  <b>
                    {samples.map((s) => (
                      <p className={styles.selection}>
                        <div className={styles.bullet}>{"\u25CF"}</div>
                        <p>{s}</p>
                      </p>
                    ))}
                  </b>
                </div>
                <div className={styles["summary-field"]}>
                  <p className={styles["selection-name"]}>Phenotypes</p>
                  <b>
                    {phenotypes?.length > 0
                      ? phenotypes.map((p) => (
                          <p className={styles.selection}>
                            <div className={styles.bullet}>{"\u25CF"}</div>
                            <p>
                              {p.ontologyId}: {p.name}
                            </p>
                          </p>
                        ))
                      : "No phenotypes selected"}
                  </b>
                </div>
                <div className={styles["summary-field"]}>
                  <p className={styles["selection-name"]}>ACMG Rules</p>
                  <b>
                    {acmgRules?.length > 0
                      ? acmgRules.map((a) => (
                          <p className={styles.selection}>
                            <div className={styles.bullet}>{"\u25CF"}</div>
                            <p>{a.name}</p>
                          </p>
                        ))
                      : "No rules selected"}
                  </b>
                </div>
              </div>
              <div className={styles.buttons}>
                <Button text="Save VUS" onClick={saveVus} />
                <Button
                  text="Modify VUS"
                  onClick={() => setIsSummaryModelOpen(false)}
                />
              </div>
            </div>
          </Modal>
        )}
      </div>
    );
  }

  function checkChrPositionValidity(val: string) {
    const isValid = val.length > 0;

    if (isValid) {
      setChromosomePositionErrorMsg("");
    } else {
      setChromosomePositionErrorMsg("Please input a chromosome position");
    }

    return isValid;
  }

  function checkAlleleValidity(val: string, isRef: boolean) {
    let isValid = /^[G|A|C|T]+$/gim.test(val);

    if (!isValid && !isRef) {
      isValid = val === "/";
    }

    if (!isValid) {
      if (isRef) {
        setRefAlleleErrorMsg(
          "Invalid allele. Must contain only 'G', 'A', 'C' and 'T'!"
        );
      } else {
        setAltAlleleErrorMsg(
          "Invalid allele. Must contain only 'G', 'A', 'C' and 'T'! If there is no reference allele input '/'."
        );
      }
    } else {
      if (isRef) {
        setRefAlleleErrorMsg("");
      } else {
        setAltAlleleErrorMsg("");
      }
    }

    return isValid;
  }

  function checkGeneValidity(val: string) {
    let isValid = val.length > 0;

    if (isValid) {
      setIsValidatingGeneInput(true);
      props.vusService?.verifyGene({ geneName: val }).then((res) => {
        if (res.isSuccess) {
          setGeneId(res.geneId);
          setGeneInputErrorMsg("");
        } else {
          setGeneInput(undefined);
          setGeneInputErrorMsg("Please enter a valid gene");
          isValid = false;
        }
        setIsValidatingGeneInput(false);
      });
    }

    return isValid;
  }

  function removeGene() {
    setGeneId(undefined);
    setGeneInput(undefined);
  }

  function removeSample(sample: string) {
    if (samples.some((s) => s === sample)) {
      const updatedSelection = samples.filter((s) => s !== sample);

      setSamples(updatedSelection);
    }
  }

  function addSample(sample: string) {
    if (sample.trim().length > 0 && !samples.some((s) => s === sample)) {
      const updatedSelection = samples.concat(sample);

      setSamples(updatedSelection);

      //clear inputted sample
      (document.getElementById("sample-input") as HTMLInputElement).value = "";
    }
  }

  function saveVus() {
    let vus: IVusUpload = {
      chromosome: chromosome,
      chromosomePosition: chromosomePosition.toString(),
      refAllele: refAllele.toUpperCase(),
      altAllele: altAllele.toUpperCase(),
      genotype:
        genotype === "Heterozygous"
          ? `${refAllele.toUpperCase()}/${altAllele.toUpperCase()}`
          : `${altAllele.toUpperCase()}/${altAllele.toUpperCase()}`,
      type: type,
      gene: geneInput.toUpperCase(),
      geneId: geneId,
      classification: classification.replace(" ", "_").toUpperCase(),
      samples: samples,
      phenotypes: phenotypes ?? [],
      acmgRules: acmgRules ?? [],
    };

    setIsLoading(true);

    props.vusService.uploadVus({ vus: vus }).then((res) => {
      if (res.isSuccess) {
        window.location.href = `/vus/${res.vusList[0].id}`;
      } else {
        //TODO: handle error
      }
    });
  }
};

export default VusUploadPage;
