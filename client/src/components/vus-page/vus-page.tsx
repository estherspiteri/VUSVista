import React, { useState } from "react";
import styles from "./vus-page.module.scss";
import { IVus } from "../../models/view-vus.model";
import VusInfo from "./vus-info/vus-info";
import { IAcmgRule } from "../../models/acmg-rule.model";
import { VusService } from "../../services/vus/vus.service";
import { Link } from "react-router-dom";
import Icon from "../../atoms/icons/icon";
import Button from "../../atoms/button/button";
import Modal from "../../atoms/modal/modal";
import Loader from "../../atoms/loader/loader";

type VusPageProps = {
  vus: IVus;
  acmgRules: IAcmgRule[];
  vusService?: VusService;
};

const VusPage: React.FunctionComponent<VusPageProps> = (
  props: VusPageProps
) => {
  const [isDeleteModalVisible, setIsDeleteModalVisible] = useState(false);
  const [isDeletingVariant, setIsDeletingVariant] = useState(false);

  return (
    <div className={styles["vus-page-container"]}>
      <div className={styles["title-container"]}>
        <div className={styles["title-section"]}>
          <div className={styles.title}>Variant Information</div>
          <div className={styles["option-container"]}>
            <div className={styles["option-btn"]}>
              <span>Variant Options</span>
              <Icon name="options" width={12} height={12} />
            </div>
            <div className={styles.options}>
              <Link to={`/review/${props.vus.id}`} className={styles.option}>
                New Classification Review
              </Link>
              <Link
                to={`/review-history/${props.vus.id}`}
                className={styles.option}
              >
                Classification Review History
              </Link>
              {props.vus.numOfPublications > 0 && (
                <Link
                  to={`/publication-view/${props.vus.id}`}
                  className={styles.option}
                >
                  View Publications
                </Link>
              )}
              <Button
                text="Delete Variant"
                icon="bin"
                className={styles["delete-btn"]}
                onClick={() => setIsDeleteModalVisible(true)}
              />
            </div>
          </div>
        </div>
        <div className={styles.description}>
          <p>Below you can find information about the selected variant.</p>
        </div>
      </div>

      <VusInfo
        vus={props.vus}
        acmgRules={props.acmgRules}
        vusService={props.vusService}
      />

      {isDeleteModalVisible && (
        <Modal
          modalContainerStyle={styles["confirm-delete-modal"]}
          isClosable={!isDeletingVariant}
          onCloseIconClickCallback={() => setIsDeleteModalVisible(false)}
        >
          <div className={styles["confirm-delete-modal-content"]}>
            <p>
              Are you sure you want to delete Variant&nbsp;
              <b>{props.vus.id}</b> ?
            </p>
            <div className={styles["option-btns"]}>
              <Button
                text="Yes, delete it"
                onClick={deleteVariant}
                className={styles["delete-btn"]}
                disabled={isDeletingVariant}
              />
              <Button
                text="No, return to variant page"
                onClick={() => setIsDeleteModalVisible(false)}
                disabled={isDeletingVariant}
              />
            </div>
          </div>
          {isDeletingVariant && <Loader />}
        </Modal>
      )}
    </div>
  );

  function deleteVariant() {
    setIsDeletingVariant(true);

    props.vusService
      .deleteVariant({ variantId: props.vus.id.toString() })
      .then((res) => {
        if (res.isSuccess) {
          window.location.href = "/view-vus";
        }

        setIsDeletingVariant(false);
      });
  }
};

export default VusPage;
