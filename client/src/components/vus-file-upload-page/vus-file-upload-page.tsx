import React, {
  useCallback,
  useContext,
  useEffect,
  useRef,
  useState,
} from "react";
import styles from "./vus-file-upload-page.module.scss";
import { VusService } from "../../services/vus/vus.service";
import { ErrorCode, FileRejection, useDropzone } from "react-dropzone";
import uploadGif from "./upload.gif";
import { Banner } from "../../atoms/banner/banner";
import Text from "../../atoms/text/text";
import Loader from "../../atoms/loader/loader";
import { INonExistingGene, IVusGene } from "../../models/vus-file-upload.model";
import Button from "../../atoms/button/button";
import Icon from "../../atoms/icons/icon";
import Modal, { ModalRef } from "../../atoms/modal/modal";
import { AppContext } from "../../app-context";
import { saveAs } from "file-saver";
type VusFileUploadPageProps = {
  vusService?: VusService;
};

type GeneNotFoundSelection = {
  geneNotFound: INonExistingGene;
  geneSelected: { name?: string; id?: number };
  errorMsg: string;
  isValidating: boolean;
};

const VusFileUploadPage: React.FunctionComponent<VusFileUploadPageProps> = (
  props: VusFileUploadPageProps
) => {
  const [file, setFile] = useState<File | undefined>(undefined);
  const [isProcessing, setIsProcessing] = useState(false);
  const [isFileProcessed, setIsFileProcessed] = useState(false);

  const [multipleGenes, setMultipleGenes] = useState<IVusGene[]>(undefined);
  const [multipleGenesSelection, setMultipleGenesSelection] = useState<
    (string | undefined)[] | undefined
  >(undefined);
  const [multipleGenesSelectionErrorMsg, setMultipleGenesSelectionErrorMsg] =
    useState("");

  const [genesNotFoundSelection, setGenesNotFoundSelection] =
    useState<GeneNotFoundSelection[]>(undefined);
  const [genesNotFoundSelectionErrorMsg, setGenesNotFoundSelectionErrorMsg] =
    useState("");

  const { taskIds, setTaskIds } = useContext(AppContext);

  const [errorMsg, setErrorMsg] = useState("");

  const multipleGenesModalRef = useRef<ModalRef>(null);
  const genesNotFoundModalRef = useRef<ModalRef>(null);

  const onFileDrop = useCallback(
    (acceptedFiles: File[], fileRejections: FileRejection[]) => {
      //ensure only one file dropped
      if (fileRejections.length > 1) {
        setErrorMsg("Too many files!");
      }
      //ensure file is a spreadsheet
      else if (
        fileRejections.length === 1 &&
        fileRejections[0].errors
          .map((e) => e.code)
          .includes(ErrorCode.FileInvalidType)
      ) {
        setErrorMsg("File must be a spreadsheet!");
      } else {
        setFile(acceptedFiles[0]);
      }
    },
    []
  );

  const { getRootProps, getInputProps, open, isDragActive } = useDropzone({
    onDrop: onFileDrop,
    maxFiles: 1,
    accept: {
      "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet": [], // XLSX
      "application/vnd.ms-excel": [], // XLS
    },
    noClick: true,
  });

  useEffect(() => {
    if (isDragActive && errorMsg.length > 0) {
      setErrorMsg("");
    }
  }, [isDragActive, errorMsg]);

  useEffect(() => {
    const areAllGenesSelected = !multipleGenesSelection?.includes(undefined);

    if (areAllGenesSelected) {
      setMultipleGenesSelectionErrorMsg("");
    }
  }, [multipleGenesSelection]);

  return (
    <div className={styles["vus-file-upload-page-container"]}>
      <div className={styles.title}>VUS File Upload</div>
      <div className={styles.description}>
        {file ? (
          <p>Process your uploaded file.</p>
        ) : (
          <>
            <p style={{ lineHeight: "28px" }}>
              Drag 'n' drop or click on the button to select your VUS file from
              the file directory. Your VUS file should be a spreadheet with the
              same headers as in this template:&nbsp;
              <span
                style={{
                  color: "#008080",
                  cursor: "pointer",
                  fontWeight: "bold",
                }}
                onClick={downloadTemplate}
              >
                template.xlsx
              </span>
              .
            </p>
          </>
        )}
      </div>
      <div className={styles["drag-and-drop-wrapper"]}>
        {!file ? (
          <div
            {...getRootProps({
              className: styles["drag-and-drop-container"],
            })}
          >
            <input {...getInputProps()} />
            <div className={styles["drag-and-drop"]}>
              <div className={styles["dotted-drag-and-drop"]}>
                {isDragActive ? (
                  <img src={uploadGif} width={200} />
                ) : (
                  <>
                    <Icon
                      name="file-upload"
                      width={110}
                      height={110}
                      stroke="#008080"
                    />
                    <p>Drag your VUS spreadsheet file here</p>
                    <p>OR</p>
                    <Button text="Browse files" onClick={open} />
                  </>
                )}
              </div>
            </div>
          </div>
        ) : (
          <div className={styles["uploaded-file-wrapper"]}>
            <div className={styles["uploaded-file-container"]}>
              <div className={styles["uploaded-file"]}>
                <span className={styles["file-name"]}>{file.name}</span>
                <div className={styles.icon} onClick={() => setFile(undefined)}>
                  {isProcessing ||
                  (errorMsg.length === 0 && multipleGenes?.length > 0) ? (
                    <Loader width={16} thickness={2} />
                  ) : isFileProcessed ? (
                    <Icon name="checkmark" fill="#008080" />
                  ) : (
                    <Icon name="close" stroke="#008080" />
                  )}
                </div>
              </div>
            </div>
            {isFileProcessed && (
              <p style={{ color: "#008080" }}>
                Your file is currently being processed. You will be notified
                when it has been uploaded.
              </p>
            )}
            <Button
              text={
                isFileProcessed
                  ? "Upload new file"
                  : isProcessing
                  ? "Processing"
                  : "Process file"
              }
              icon="file-process"
              disabled={isProcessing}
              onClick={() =>
                isFileProcessed ? newFileUpload() : processFile()
              }
            />
          </div>
        )}
      </div>
      {errorMsg.length > 0 && (
        <Banner isClosable={true}>
          <p className={styles.errorMsg}>{errorMsg}</p>
        </Banner>
      )}
      {!file && (
        <div className={styles.important}>
          <span className={styles["important-title"]}>IMPORTANT</span>
          <ul>
            <li>
              <span>
                <b>Sample Ids</b>: need to be separated by a <b>comma</b>
              </span>
            </li>
            <li>
              <span>
                <b>Gene</b>: if a variant has multiple genes they need to be
                separated by a <b>comma</b>
              </span>
            </li>
            <li>
              <span>
                <b>ACMG rules</b>: need to be separated by a <b>comma</b>
              </span>
            </li>
            <li>
              <span>
                <b>Literature links</b>: need to be separated by a&nbsp;
                <b>pipeline</b>, i.e. |
              </span>
            </li>
            <li>
              <span>
                RSID header needs to end with an <b>underscore</b>, i.e. _
              </span>
            </li>
          </ul>
        </div>
      )}

      {multipleGenes && multipleGenes.length > 0 && (
        <Modal
          ref={multipleGenesModalRef}
          title={`Please select one Gene for each of the following Variants:`}
          isClosable={true}
          onCloseIconClickCallback={newFileUpload}
        >
          <div className={styles["modal-content-wrapper"]}>
            <div className={styles["selection-content"]}>
              <div className={styles["variant-info"]}>
                <div className={styles.header}>
                  <div className={styles.field}>Locus</div>
                  <div className={styles.field}>Type</div>
                  <div className={styles.field}>Reference</div>
                  <div className={styles.field}>Alternate Allele</div>
                  <div
                    className={`${styles.field} ${styles["multiple-genes"]}`}
                  >
                    Gene Selection
                  </div>
                </div>
                <div className={styles.content}>
                  {multipleGenes.map((x, index) => {
                    return (
                      <div className={styles.info}>
                        <div className={styles.field}>{x.vus.locus}</div>
                        <div className={styles.field}>{x.vus.type}</div>
                        <div className={styles.field}>{x.vus.refAllele}</div>
                        <div className={styles.field}>{x.vus.altAllele}</div>
                        <div
                          className={`${styles.field} ${styles["multiple-genes"]}`}
                        >
                          {x.genes.map((gene) => {
                            return (
                              <label className={styles["gene-label"]}>
                                <input
                                  type="radio"
                                  name={`gene-selection-${index}`}
                                  value={gene}
                                  checked={
                                    multipleGenesSelection[index] === gene
                                  }
                                  onChange={(e) =>
                                    setMultipleGenesSelection(
                                      multipleGenesSelection.map(
                                        (selection, i) => {
                                          if (i === index) {
                                            return gene;
                                          } else return selection;
                                        }
                                      )
                                    )
                                  }
                                  className={styles.radio}
                                />
                                {gene}
                              </label>
                            );
                          })}
                        </div>
                      </div>
                    );
                  })}
                </div>
              </div>
            </div>
            {multipleGenesSelectionErrorMsg?.length > 0 && (
              <span className={styles["error-msg"]}>
                {multipleGenesSelectionErrorMsg}
              </span>
            )}
            <Button
              text="Continue processing file"
              icon="file-process"
              onClick={saveMultipleGenesSelection}
            />
          </div>
        </Modal>
      )}

      {genesNotFoundSelection && genesNotFoundSelection.length > 0 && (
        <Modal
          ref={genesNotFoundModalRef}
          title={`Please input a valid Gene for each of the following Variants:`}
          isClosable={true}
          onCloseIconClickCallback={newFileUpload}
        >
          <p style={{ marginBottom: "16px" }}>
            The genes for the following variants were not found in our database.
            Kindly input a valid gene and click on the search icon for it to be
            validated. Once a gene is validated, the text box will no longer be
            visible. Check out gene aliases at&nbsp;
            <a
              href="https://www.genecards.org/"
              target="_blank"
              rel="noopener noreferrer"
              color="#008080"
            >
              Gene Cards
            </a>
          </p>
          <div className={styles["modal-content-wrapper"]}>
            <div className={styles["selection-content"]}>
              <div className={styles["not-found-selected-gene-info"]}>
                <div className={styles.header}>
                  <div className={styles.field}>Locus</div>
                  <div className={styles.field}>Type</div>
                  <div className={styles.field}>Reference</div>
                  <div className={styles.field}>Alternate Allele</div>
                  <div className={styles.field}>Invalid Gene</div>
                  <div
                    className={`${styles.field} ${styles["not-found-selected-gene"]}`}
                  >
                    Valid Gene Check
                  </div>
                </div>
                <div className={styles.content}>
                  {genesNotFoundSelection.map((x, genesNotFoundIndex) => {
                    return (
                      <div className={styles.info}>
                        <div className={styles.field}>
                          {x.geneNotFound.vus.locus}
                        </div>
                        <div className={styles.field}>
                          {x.geneNotFound.vus.type}
                        </div>
                        <div className={styles.field}>
                          {x.geneNotFound.vus.refAllele}
                        </div>
                        <div className={styles.field}>
                          {x.geneNotFound.vus.altAllele}
                        </div>
                        <div className={styles.field}>
                          {x.geneNotFound.gene}
                        </div>
                        <div
                          className={`${styles.field} ${styles["not-found-selected-gene"]}`}
                        >
                          {genesNotFoundSelection[genesNotFoundIndex]
                            .geneSelected.id ? (
                            <div className={styles["gene-selected"]}>
                              <span>
                                {genesNotFoundSelection[
                                  genesNotFoundIndex
                                ].geneSelected.name.toUpperCase()}
                              </span>
                              <Icon
                                name="close"
                                onClick={() => removeGene(genesNotFoundIndex)}
                                stroke="#008080"
                                cursor="pointer"
                              />
                            </div>
                          ) : (
                            <div className={styles["gene-input"]}>
                              <Text
                                value={
                                  genesNotFoundSelection[genesNotFoundIndex]
                                    .geneSelected.name
                                }
                                disabled={
                                  genesNotFoundSelection[genesNotFoundIndex]
                                    .isValidating
                                }
                                onChange={(e) =>
                                  setGenesNotFoundSelection(
                                    genesNotFoundSelection.map(
                                      (selection, i) => {
                                        if (i === genesNotFoundIndex) {
                                          return {
                                            ...genesNotFoundSelection[i],
                                            errorMsg: "",
                                            geneSelected: {
                                              name: e.currentTarget.value,
                                            },
                                          };
                                        } else return selection;
                                      }
                                    )
                                  )
                                }
                                errorMsg={
                                  genesNotFoundSelection[genesNotFoundIndex]
                                    .errorMsg
                                }
                              />
                              <Icon
                                name="search"
                                onClick={() => {
                                  if (
                                    !genesNotFoundSelection[genesNotFoundIndex]
                                      .isValidating
                                  ) {
                                    checkGeneValidity(genesNotFoundIndex);
                                  }
                                }}
                                stroke="#008080"
                                fill="#fff"
                                cursor="pointer"
                              />
                            </div>
                          )}
                        </div>
                      </div>
                    );
                  })}
                </div>
              </div>
            </div>
            {genesNotFoundSelectionErrorMsg?.length > 0 && (
              <span className={styles["error-msg"]}>
                {genesNotFoundSelectionErrorMsg}
              </span>
            )}
            <Button
              text="Continue processing file"
              icon="file-process"
              onClick={saveGenesNotFoundSelection}
            />
          </div>
        </Modal>
      )}
    </div>
  );

  function saveMultipleGenesSelection() {
    const areAllGenesSelected = !multipleGenesSelection.includes(undefined);

    if (areAllGenesSelected) {
      checkGenesExist();
      multipleGenesModalRef.current.closeModal();
    } else {
      setMultipleGenesSelectionErrorMsg(
        "Please select a single gene for each of the above variants!"
      );
    }
  }

  function saveGenesNotFoundSelection() {
    const areThereUnselectedGenes = genesNotFoundSelection.find(
      (g) => !g.geneSelected.id
    );

    if (areThereUnselectedGenes) {
      setGenesNotFoundSelectionErrorMsg(
        "Please input and verfiy genes for each of the above variants!"
      );
    } else {
      storeAndVerifyFile();
      genesNotFoundModalRef.current.closeModal();
    }
  }

  function newFileUpload() {
    setFile(undefined);

    setIsProcessing(false);
    setIsFileProcessed(false);
    setMultipleGenes(undefined);
    setMultipleGenesSelection(undefined);
    setMultipleGenesSelectionErrorMsg("");
    setGenesNotFoundSelection(undefined);
    setGenesNotFoundSelectionErrorMsg("");

    setErrorMsg("");
  }

  function storeAndVerifyFile() {
    props.vusService
      ?.storeAndVerifyVusFile({
        vusFile: file,
        genesNotFoundSelection: genesNotFoundSelection?.map((g) => {
          return { index: g.geneNotFound.index, gene: g.geneSelected.name };
        }),
        multipleGenesSelection: multipleGenesSelection?.map((g, i) => {
          return { index: multipleGenes[i].index, gene: g };
        }),
      })
      .then((res) => {
        setIsFileProcessed(true);
        setIsProcessing(false);

        if (!res.isSuccess) {
          setErrorMsg("Failed to succesfully process file. Please try again!");
        } else {
          setMultipleGenes(undefined);
          setMultipleGenesSelection(undefined);
          setTaskIds(taskIds.concat(res.taskId));
        }
      });
  }

  function checkGenesExist() {
    // check if all genes can be found in DB
    props.vusService
      .checkFileForValidGenes({
        vusFile: file,
        multipleGenesSelection: multipleGenesSelection?.map((g, i) => {
          return { index: multipleGenes[i].index, gene: g };
        }),
      })
      .then((res) => {
        if (res.genesNotInDb && res.genesNotInDb.length > 0) {
          setIsFileProcessed(true);

          // initialise gene selection array
          setGenesNotFoundSelection(
            res.genesNotInDb.map((g) => {
              return {
                geneNotFound: g,
                geneSelected: { name: undefined, id: undefined },
                errorMsg: "",
                isValidating: false,
              };
            })
          );
          return;
        } else {
          storeAndVerifyFile();
        }
      });
  }

  function processFile() {
    setIsProcessing(true);
    setErrorMsg("");

    props.vusService
      .checkFileForMultipleGenes({
        vusFile: file,
      })
      .then((res) => {
        if (res.multipleGenes && res.multipleGenes.length > 0) {
          setIsFileProcessed(true);

          // populate gene selection array with undefined
          setMultipleGenesSelection(
            Array.apply(undefined, Array(res.multipleGenes.length))
          );

          setMultipleGenes(res.multipleGenes);
          return;
        } else {
          checkGenesExist();
        }
      });
  }

  function checkGeneValidity(index: number) {
    let isValid = false;

    // only validate if the gene has no errors
    if (genesNotFoundSelection[index].errorMsg.length === 0) {
      const val = genesNotFoundSelection[index].geneSelected.name;

      isValid = val?.length > 0;

      if (isValid) {
        //set is validating to true
        setGenesNotFoundSelection(
          genesNotFoundSelection.map((selection, i) => {
            if (i === index) {
              return {
                ...genesNotFoundSelection[i],
                isValidating: true,
              };
            } else return selection;
          })
        );

        // verify the gene
        props.vusService?.verifyGene({ geneName: val }).then((res) => {
          setGenesNotFoundSelection(
            genesNotFoundSelection.map((selection, i) => {
              if (i === index) {
                if (res.isSuccess) {
                  return {
                    ...genesNotFoundSelection[i],
                    geneSelected: { name: val, id: res.geneId },
                    errorMsg: "",
                    isValidating: false,
                  };
                } else {
                  return {
                    ...genesNotFoundSelection[i],
                    geneSelected: { name: undefined },
                    errorMsg: "Please enter a valid gene",
                    isValidating: false,
                  };
                }
              } else return selection;
            })
          );
        });
      }
    }

    return isValid;
  }

  function removeGene(index: number) {
    setGenesNotFoundSelection(
      genesNotFoundSelection.map((selection, i) => {
        if (i === index) {
          return {
            ...genesNotFoundSelection[i],
            geneSelected: { id: undefined, name: undefined },
          };
        } else return selection;
      })
    );
  }

  function downloadTemplate() {
    const fileUrl = process.env.PUBLIC_URL + "/assets/template.xlsx";

    // Fetch the file
    fetch(fileUrl)
      .then((response) => {
        if (!response.ok) {
          throw new Error("Network response was not ok");
        }
        return response.blob();
      })
      .then((blob) => {
        // Use file-saver to save the file
        saveAs(blob, "template.xlsx");
      })
      .catch((error) => {
        console.error("There was a problem with the fetch operation:", error);
      });
  }
};

export default VusFileUploadPage;
