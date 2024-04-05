import React, { useState } from "react";
import styles from "./vus-upload-page.module.scss";
import Text from "./../../atoms/text/text";
import VusUploadField from "./vus-upload-field/vus-upload-field";
import Button from "../../atoms/button/button";
import SamplePhenotypeSelection from "../sample-page/sample-phenotype-selection/sample-phenotype-selection";
import { IHPOTerm } from "../../services/sample/sample.dto";
import { SampleService } from "../../services/sample/sample.service";
import Icon from "../../atoms/icons/icon";
import { IGene } from "../../models/gene.model";
import Dropdown from "../../atoms/dropdown/dropdown";
import { VusService } from "../../services/vus/vus.service";
import Loader from "../../atoms/loader/loader";
import Modal from "../../atoms/modal/modal";

type VusUploadPageProps = {
  sampleService?: SampleService;
  vusService?: VusService;
};

//TODO: In the end show variant summary - let user confirm on cancel submittion to continue editing
const VusUploadPage: React.FunctionComponent<VusUploadPageProps> = (
  props: VusUploadPageProps
) => {
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

  const [geneInput, setGeneInput] = useState<string | undefined>(undefined);
  const [isValidatingGeneInput, setIsValidatingGeneInput] = useState(false);
  const [geneInputErrorMsg, setGeneInputErrorMsg] = useState("");
  const [geneId, setGeneId] = useState<number | undefined>(undefined);

  const [classification, setClassification] = useState<string | undefined>(
    undefined
  );
  const [classificationErrorMsg, setClassificationErrorMsg] = useState("");

  const [samples, setSamples] = useState<string[]>([]);
  const [samplesErrorMsg, setSamplesErrorMsg] = useState("");

  const [phenotypes, setPhenotypes] = useState<IHPOTerm[]>(undefined);

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

  const isGeneValid = geneId !== undefined;

  const isClassificationValid =
    classification !== undefined && classificationErrorMsg.length === 0;

  const areSamplesValid = samples.length > 0;

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
                  <span>{geneInput}</span>
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

        {/** Samples */}
        <VusUploadField title="Samples" showCheckMark={areSamplesValid}>
          <div className={styles["field-content-container"]}>
            <div className={styles["field-content"]}>
              <span>
                List the samples that have this variant below. Press 'Enter' to
                add the sample to the list.
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
                          <span className={styles["selected-sample"]}>{s}</span>
                        </p>
                      );
                    })}
                </div>
              </div>
            </div>
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
      </div>

      <Button
        className={styles.btn}
        text="Save variant"
        disabled={
          !isChrValid ||
          !areAllelesValid ||
          !isGeneValid ||
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
                <p className={styles.selection}>
                  <div className={styles.bullet}>{"\u25CF"}</div>
                  <b>
                    chr{chromosome}:{chromosomePosition}
                  </b>
                </p>
              </div>
              <div className={styles["summary-field"]}>
                <p className={styles["selection-name"]}>Reference allele</p>
                <p className={styles.selection}>
                  <div className={styles.bullet}>{"\u25CF"}</div>
                  <b>{refAllele.toUpperCase()}</b>
                </p>
              </div>
              <div className={styles["summary-field"]}>
                <p className={styles["selection-name"]}>Alternate allele</p>
                <p className={styles.selection}>
                  <div className={styles.bullet}>{"\u25CF"}</div>
                  <b>{altAllele.toUpperCase()}</b>
                </p>
              </div>
              <div className={styles["summary-field"]}>
                <p className={styles["selection-name"]}>Gene</p>
                <p className={styles.selection}>
                  <div className={styles.bullet}>{"\u25CF"}</div>
                  <b>{geneInput}</b>
                </p>
              </div>
              <div className={styles["summary-field"]}>
                <p className={styles["selection-name"]}>Classification</p>
                <p className={styles.selection}>
                  <div className={styles.bullet}>{"\u25CF"}</div>
                  <b>{classification}</b>
                </p>
              </div>

              <div className={styles["summary-field"]}>
                <p className={styles["selection-name"]}>Sample Ids</p>
                <p className={styles.selection}>
                  <b>
                    {samples.map((s) => (
                      <div>
                        <div className={styles.bullet}>{"\u25CF"}</div>
                        <p>{s}</p>
                      </div>
                    ))}
                  </b>
                </p>
              </div>
              <div className={styles["summary-field"]}>
                <p className={styles["selection-name"]}>Phenotypes</p>
                <div className={styles.selection}>
                  <b>
                    {phenotypes?.length > 0
                      ? phenotypes.map((p) => (
                          <div>
                            <div className={styles.bullet}>{"\u25CF"}</div>
                            <p>
                              {p.ontologyId}: {p.name}
                            </p>
                          </div>
                        ))
                      : "No phenotypes selected"}
                  </b>
                </div>
              </div>
            </div>
            <div className={styles.buttons}>
              <Button text="Save VUS" />
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
    if (!samples.some((s) => s === sample)) {
      const updatedSelection = samples.concat(sample);

      setSamples(updatedSelection);

      //clear inputted sample
      (document.getElementById("sample-input") as HTMLInputElement).value = "";
    }
  }
};

export default VusUploadPage;
