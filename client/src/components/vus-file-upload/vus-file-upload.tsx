import React, { useRef, useState } from "react";
import styles from "./vus-file-upload.module.scss";
import { VusService } from "../../services/vus/vus.service";

type VusFileUploadProps = {
  vusService?: VusService;
};

const VusFileUpload: React.FunctionComponent<VusFileUploadProps> = (
  props: VusFileUploadProps
) => {
  const [isUploadSuccess, setIsUploadSuccess] = useState(undefined);
  const inputRef = useRef<HTMLInputElement>(null);

  return (
    <div className={styles["vus-file-upload-container"]}>
      <input ref={inputRef} type="file" />
      <button onClick={handleFileUpload}>Load VUS file</button>
      <h1>Is File Upload Success: {isUploadSuccess}</h1>
    </div>
  );

  function handleFileUpload(e: any) {
    e.preventDefault();

    props.vusService
      ?.storeAndVerifyVusFile({
        vusFile: inputRef.current.files[0],
      })
      .then((res) => {
        setIsUploadSuccess(res.isSuccess);
      });
  }
};

export default VusFileUpload;
