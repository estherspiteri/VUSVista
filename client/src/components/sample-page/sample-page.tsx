import React, { useState } from "react";
import styles from "./sample-page.module.scss";
import { ISample } from "../../models/view-samples.model";
import SampleInfo from "./sample-info/sample-info";
import { SampleService } from "../../services/sample/sample.service";
import Button from "../../atoms/button/button";
import Modal from "../../atoms/modal/modal";
import Loader from "../../atoms/loader/loader";

type SamplePageProps = {
  sample: ISample;
  sampleService: SampleService;
};

const SamplePage: React.FunctionComponent<SamplePageProps> = (
  props: SamplePageProps
) => {
  const [isDeleteModalVisible, setIsDeleteModalVisible] = useState(false);
  const [isDeletingSample, setIsDeletingSample] = useState(false);

  return (
    <div className={styles["sample-page-container"]}>
      <div className={styles["title-container"]}>
        <div className={styles.title}>Sample Information</div>
        <Button
          text="Delete Sample"
          icon="bin"
          className={styles["delete-btn"]}
          onClick={() => setIsDeleteModalVisible(true)}
        />
      </div>
      <div className={styles.description}>
        <p>Below you can find information about the selected sample.</p>
      </div>

      <SampleInfo sample={props.sample} sampleService={props.sampleService} />

      {isDeleteModalVisible && (
        <Modal
          modalContainerStyle={styles["confirm-delete-modal"]}
          isClosable={!isDeletingSample}
          onCloseIconClickCallback={() => setIsDeleteModalVisible(false)}
        >
          <div className={styles["confirm-delete-modal-content"]}>
            <p>
              Are you sure you want to delete Sample&nbsp;
              <b>{props.sample.sampleId}</b> ?
            </p>
            <div className={styles["option-btns"]}>
              <Button
                text="Yes, delete it"
                onClick={deleteSample}
                className={styles["delete-btn"]}
                disabled={isDeletingSample}
              />
              <Button
                text="No, return to sample page"
                onClick={() => setIsDeleteModalVisible(false)}
                disabled={isDeletingSample}
              />
            </div>
          </div>
          {isDeletingSample && <Loader />}
        </Modal>
      )}
    </div>
  );

  function deleteSample() {
    setIsDeletingSample(true);

    props.sampleService
      .deleteSample({ sampleId: props.sample.sampleId })
      .then((res) => {
        if (res.isSuccess) {
          window.location.href = "/view-samples";
        }

        setIsDeletingSample(false);
      });
  }
};

export default SamplePage;
