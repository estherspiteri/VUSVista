import React, { useCallback, useEffect, useState } from "react";
import styles from "./vus-file-upload-page.module.scss";
import { VusService } from "../../services/vus/vus.service";
import { ErrorCode, FileRejection, useDropzone } from "react-dropzone";
import uploadGif from "./upload.gif";
import { Banner } from "../banner/banner";
import Loader from "../loader/loader";
import { IVus } from "../../models/view-vus.model";
import VusTable from "../view-vus-page/vus-table/vus-table";

type VusFileUploadPageProps = {
  vusService?: VusService;
};

const VusFileUploadPage: React.FunctionComponent<VusFileUploadPageProps> = (
  props: VusFileUploadPageProps
) => {
  const [areRsidsRetrieved, setAreRsidsRetrieved] = useState(false);
  const [isClinvarAccessed, setIsClinvarAccessed] = useState(false);

  const [file, setFile] = useState<File | undefined>(undefined);
  const [isProcessing, setIsUploading] = useState(false);
  const [isFileProcessed, setIsFileUploaded] = useState(false);
  const [vusList, setVusList] = useState<IVus[]>(undefined);

  const [errorMsg, setErrorMsg] = useState("");

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
                    <svg
                      width="198"
                      height="198"
                      viewBox="0 0 24 24"
                      fill="none"
                      xmlns="http://www.w3.org/2000/svg"
                    >
                      <path
                        d="M12 16v-6m0 0-3 2m3-2 3 2M3 6v10.8c0 1.12 0 1.68.218 2.108a2 2 0 0 0 .874.874c.427.218.987.218 2.105.218h11.606c1.118 0 1.677 0 2.104-.218.377-.192.683-.498.875-.874C21 18.48 21 17.92 21 16.8V9.2c0-1.12 0-1.68-.218-2.108a2 2 0 0 0-.874-.874C19.48 6 18.92 6 17.8 6H12M3 6h9M3 6a2 2 0 0 1 2-2h3.675c.489 0 .734 0 .964.055.204.05.399.13.578.24.202.124.375.297.72.643L12 6"
                        stroke="teal"
                        stroke-width="2"
                        stroke-linecap="round"
                        stroke-linejoin="round"
                      />
                    </svg>

                    <p>Drag your VUS spreadsheet file here</p>
                    <p>OR</p>
                    <button onClick={open} className={styles.btn}>
                      Browse files
                    </button>
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
                  {isProcessing ? (
                    <Loader width={16} thickness={2} />
                  ) : isFileProcessed ? (
                    <svg
                      viewBox="0 0 24 24"
                      fill="none"
                      xmlns="http://www.w3.org/2000/svg"
                    >
                      <path
                        fill-rule="evenodd"
                        clip-rule="evenodd"
                        d="M19.707 6.293a1 1 0 0 1 0 1.414L10.414 17a2 2 0 0 1-2.828 0l-4.293-4.293a1 1 0 1 1 1.414-1.414L9 15.586l9.293-9.293a1 1 0 0 1 1.414 0Z"
                        fill="#008080"
                      />
                    </svg>
                  ) : (
                    <svg
                      viewBox="0 0 24 24"
                      fill="none"
                      xmlns="http://www.w3.org/2000/svg"
                    >
                      <path
                        d="m16 16-4-4m0 0L8 8m4 4 4-4m-4 4-4 4"
                        stroke="#008080"
                        stroke-width="2"
                        stroke-linecap="round"
                        stroke-linejoin="round"
                      />
                    </svg>
                  )}
                </div>
              </div>
            </div>
            <button
              onClick={() =>
                isFileProcessed ? newFileUpload() : processFile()
              }
              className={`${styles.btn} ${isProcessing ? styles.disabled : ""}`}
            >
              <div className={styles["btn-content"]}>
                <span>
                  {isFileProcessed
                    ? "Upload new file"
                    : isProcessing
                    ? "Processing"
                    : "Process file"}
                </span>
                <div className={styles.icon}>
                  <svg viewBox="0 0 16 16" xmlns="http://www.w3.org/2000/svg">
                    <g fill="#fff">
                      <path d="M12 0H5v6h.7l.2.7.1.1V1h5v4h4v9H9l.3.5-.5.5H16V4l-4-4zm0 4V1l3 3h-3zm-6.5 7.5a1 1 0 1 1-2 0 1 1 0 0 1 2 0z" />
                      <path d="M7.9 12.4 9 12v-1l-1.1-.4c-.1-.3-.2-.6-.4-.9l.5-1-.7-.7-1 .5c-.3-.2-.6-.3-.9-.4L5 7H4l-.4 1.1c-.3.1-.6.2-.9.4l-1-.5-.7.7.5 1.1c-.2.3-.3.6-.4.9L0 11v1l1.1.4c.1.3.2.6.4.9l-.5 1 .7.7 1.1-.5c.3.2.6.3.9.4L4 16h1l.4-1.1c.3-.1.6-.2.9-.4l1 .5.7-.7-.5-1.1c.2-.2.3-.5.4-.8zm-3.4 1.1c-1.1 0-2-.9-2-2s.9-2 2-2 2 .9 2 2-.9 2-2 2z" />
                    </g>
                  </svg>
                </div>
              </div>
            </button>
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

                <VusTable vusList={vusList} showGenotype={true} />
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
    </div>
  );

  function newFileUpload() {
    setAreRsidsRetrieved(false);
    setIsClinvarAccessed(false);

    setFile(undefined);

    setIsUploading(false);
    setIsFileUploaded(false);
    setVusList(undefined);

    setErrorMsg("");
  }

  function processFile() {
    setIsUploading(true);
    setErrorMsg("");

    props.vusService
      ?.storeAndVerifyVusFile({
        vusFile: file,
      })
      .then((res) => {
        setIsUploading(false);
        setIsFileUploaded(res.isSuccess);

        if (!res.isSuccess) {
          setErrorMsg("Failed to succesfully process file. Please try again!");
        } else {
          setVusList(res.vusList);
        }

        setAreRsidsRetrieved(res.areRsidsRetrieved);
        setIsClinvarAccessed(res.isClinvarAccessed);
      });
  }
};

export default VusFileUploadPage;
