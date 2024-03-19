import React, { useCallback, useEffect, useRef, useState } from "react";
import styles from "./vus-file-upload-page.module.scss";
import { VusService } from "../../services/vus/vus.service";
import { ErrorCode, FileRejection, useDropzone } from "react-dropzone";
import uploadGif from "./upload.gif";
import { Banner } from "../../atoms/banner/banner";
import Loader from "../../atoms/loader/loader";
import { IVus } from "../../models/view-vus.model";
import VusTable from "../view-vus-page/vus-table/vus-table";
import { IVusGene } from "../../models/vus-file-upload.model";
import Button from "../../atoms/button/button";
import Icon from "../../atoms/icon/icon";
import Modal, { ModalRef } from "../../atoms/modal/modal";

type VusFileUploadPageProps = {
  vusService?: VusService;
};

const VusFileUploadPage: React.FunctionComponent<VusFileUploadPageProps> = (
  props: VusFileUploadPageProps
) => {
  // const [areRsidsRetrieved, setAreRsidsRetrieved] = useState(false);
  // const [isClinvarAccessed, setIsClinvarAccessed] = useState(false);

  const [file, setFile] = useState<File | undefined>(undefined);
  const [isProcessing, setIsProcessing] = useState(false);
  const [isFileProcessed, setIsFileUploaded] = useState(false);
  const [vusList, setVusList] = useState<IVus[]>(undefined);

  const [multipleGenes, setMultipleGenes] = useState<IVusGene[]>(undefined);
  const [multipleGenesSelection, setMultipleGenesSelection] = useState<
    (string | undefined)[] | undefined
  >(undefined);
  const [multipleGenesSelectionErrorMsg, setMultipleGenesSelectionErrorMsg] =
    useState("");

  const [errorMsg, setErrorMsg] = useState("");

  const modalRef = useRef<ModalRef>(null);

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
          <p>
            Drag 'n' drop or click on the button to select your VUS file from
            the file directory.
          </p>
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
            {isFileProcessed && vusList && (
              <div className={styles["file-content"]}>
                <div className={styles["file-summary"]}>
                  <div className={styles.summary}>
                    <p>RSIDs</p>
                    <p>
                      {
                        vusList.filter(
                          (vus) =>
                            vus.rsid !== "NORSID" && vus.rsidDbsnpVerified
                        ).length
                      }
                    </p>
                  </div>
                  <div className={`${styles.summary} ${styles["num-vus"]}`}>
                    <p>VUS</p> <p>{vusList.length}</p>
                  </div>
                  <div className={styles.summary}>
                    <p>ClinVar</p>
                    <p>
                      {
                        vusList.filter(
                          (vus) => vus.clinvarErrorMsg.length === 0
                        ).length
                      }
                    </p>
                  </div>
                </div>

                <VusTable
                  vusList={vusList}
                  showGenotype={true}
                  showZygosity={false}
                />
              </div>
            )}
          </div>
        )}
      </div>
      {errorMsg.length > 0 && (
        <Banner isClosable={true}>
          <p className={styles.errorMsg}>{errorMsg}</p>
        </Banner>
      )}
      {multipleGenes && multipleGenes.length > 0 && (
        <Modal
          ref={modalRef}
          title={`Please select one Gene for each of the following Variants:`}
        >
          <div className={styles["modal-content-wrapper"]}>
            <div className={styles["selection-content"]}>
              <div className={styles["variant-info"]}>
                <div className={styles.header}>
                  <div className={styles.field}>Locus</div>
                  <div className={styles.field}>Type</div>
                  <div className={styles.field}>Genotype</div>
                  <div className={styles.field}>Reference</div>
                  <div className={styles.field}>Alternate Allele</div>
                </div>
                <div className={styles.content}>
                  {multipleGenes.map((x) => {
                    return (
                      <div className={styles.info}>
                        <div className={styles.field}>{x.vus.locus}</div>
                        <div className={styles.field}>{x.vus.type}</div>
                        <div className={styles.field}>{x.vus.genotype}</div>
                        <div className={styles.field}>{x.vus.refAllele}</div>
                        <div className={styles.field}>{x.vus.altAllele}</div>
                      </div>
                    );
                  })}
                </div>
              </div>
              <div className={styles["gene-selection"]}>
                {multipleGenes.map((x, multipleGenesIndex) => {
                  return (
                    <div className={`${styles.field} ${styles.genes}`}>
                      {x.genes.map((gene, geneIndex) => {
                        return (
                          <label>
                            <input
                              type="radio"
                              name={`gene-selection-${multipleGenesIndex}`}
                              value={gene}
                              checked={
                                multipleGenesSelection[multipleGenesIndex] ===
                                gene
                              }
                              onChange={(e) =>
                                setMultipleGenesSelection(
                                  multipleGenesSelection.map((selection, i) => {
                                    if (i === multipleGenesIndex) {
                                      return gene;
                                    } else return selection;
                                  })
                                )
                              }
                              className={styles.radio}
                            />
                            {gene}
                          </label>
                        );
                      })}
                    </div>
                  );
                })}
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
    </div>
  );

  function saveMultipleGenesSelection() {
    const areAllGenesSelected = !multipleGenesSelection.includes(undefined);

    if (areAllGenesSelected) {
      processFile();
      modalRef.current.closeModal();
    } else {
      setMultipleGenesSelectionErrorMsg(
        "Please select a single gene for each of the above variants!"
      );
    }
  }

  function newFileUpload() {
    // setAreRsidsRetrieved(false);
    // setIsClinvarAccessed(false);

    setFile(undefined);

    setIsProcessing(false);
    setIsFileUploaded(false);
    setVusList(undefined);

    setErrorMsg("");
  }

  function processFile() {
    setIsProcessing(true);
    setErrorMsg("");

    props.vusService
      ?.storeAndVerifyVusFile({
        vusFile: file,
        multipleGenesSelection: multipleGenesSelection?.map((g, i) => {
          return { index: multipleGenes[i].index, gene: g };
        }),
      })
      .then((res) => {
        setIsFileUploaded(res.isSuccess);

        if (!res.isSuccess) {
          setIsProcessing(false);
          setErrorMsg("Failed to succesfully process file. Please try again!");
        } else {
          if (res.multipleGenes && res.multipleGenes.length > 0) {
            // populate gene selection array with undefined
            setMultipleGenesSelection(
              Array.apply(undefined, Array(res.multipleGenes.length))
            );

            setMultipleGenes(res.multipleGenes);
          } else {
            setIsProcessing(false);
            setMultipleGenes(undefined);
            setMultipleGenesSelection(undefined);
            setVusList(res.vusList);
          }
        }

        // setAreRsidsRetrieved(res.areRsidsRetrieved);
        // setIsClinvarAccessed(res.isClinvarAccessed);
      });
  }
};

export default VusFileUploadPage;
