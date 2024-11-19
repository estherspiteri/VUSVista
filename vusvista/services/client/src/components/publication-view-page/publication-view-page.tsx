import React, { useRef, useState } from "react";
import styles from "./publication-view-page.module.scss";
import PublicationPreview from "./publication-preview/publication-preview";
import { IPublicationPreview } from "../../models/publication-view.model";
import VariantSummary from "../shared/variant-summary/variant-summary";
import { IVUSSummary } from "../../models/vus-summary.model";
import DebouncedInput from "../shared/table/debounced-input/debounced-input";
import Modal from "../../atoms/modal/modal";
import CalendarDisplay from "../../atoms/calendar-display/calendar-display";
import Loader from "../../atoms/loader/loader";
import { IVariantPublicationUpdates } from "../../models/variant-publication-updates";
import { VusService, vusService } from "../../services/vus/vus.service";
import { openInNewWindow } from "../../helpers/open-links";
import Button from "../../atoms/button/button";
import Text from "../../atoms/text/text";
import Icon from "../../atoms/icons/icon";
import { convertPubDates } from "../../helpers/date-helper";
import { useReactToPrint } from "react-to-print";
import AddPublications from "../shared/add-publications/add-publications";

type PublicationViewPageProps = {
  variantId: string;
  variant: IVUSSummary;
  isPhenotypePublicationPage?: boolean;
  phenotype?: string;
  publications?: IPublicationPreview[];
  vusService?: VusService;
};

const PublicationViewPage: React.FunctionComponent<PublicationViewPageProps> = (
  props: PublicationViewPageProps
) => {
  const [filterValue, setFilterValue] = useState("");

  const [publications, setPublications] = useState(props.publications ?? []);

  const [showPublicationArchiveModal, setShowPublicationArchiveModal] =
    useState(false);
  const [publicationUpdates, setPublicationUpdates] = useState<
    IVariantPublicationUpdates[]
  >([]);
  const [datesWithUpdates, setDatesWithUpdates] = useState([]);

  const [showAddPublicationModal, setShowAddPublicationModal] = useState(false);
  const [addedUrls, setAddedUrls] = useState<string[]>([]);
  const [isAddingPublications, setIsAddingPublications] = useState(false);

  const [publicationToBeRemoved, setPublicationToBeRemoved] =
    useState<IPublicationPreview>(undefined);
  const [isRemovingPublication, setIsRemovingPublication] = useState(false);

  const printRef = useRef<HTMLDivElement>(null);

  const handlePrint = useReactToPrint({
    content: () => printRef.current,
  });

  return (
    <div className={styles["publication-view-container"]} ref={printRef}>
      <div className={styles["title-container"]}>
        <div className={styles["title-section"]}>
          <div className={styles.title}>Publications</div>
          <div className={styles["title-btns"]}>
            {!props.isPhenotypePublicationPage && (
              <Button
                text="Add Publications"
                icon="add"
                className={styles["add-publications-header-btn"]}
                onClick={() => setShowAddPublicationModal(true)}
              />
            )}
            {publications.length > 0 && (
              <Button
                text={"Save as PDF"}
                onClick={handlePrint}
                icon="save"
                className={styles["download-btn"]}
              />
            )}
          </div>
        </div>

        <div className={styles.description}>
          {props.isPhenotypePublicationPage ? (
            <>
              <p>
                Below you can find the publications for&nbsp;
                <b>VUS with Id {props.variantId}</b> in relation to the
                phenotype&nbsp;
                <b>{props.phenotype}</b>.
              </p>
              <p className={styles.instructions}>
                Click on a publication title to view a summary of the respective
                publication. You can view the publications in a new window by
                clicking on the button found on the right-side of each title.
              </p>
            </>
          ) : (
            <>
              <p>
                Below you can find the publications for VUS with&nbsp;
                <b>Id {props.variantId}</b>. Publications with the head icon
                next to their title were inputted manually by a scientific
                member of staff.
              </p>
              <p className={styles.instructions}>
                Click on a publication title to view a summary of the respective
                publication. You can view the publications in a new window by
                clicking on the button found on the right-side of each title.
              </p>
            </>
          )}
        </div>
      </div>

      {publications && (
        <div className={styles["publications-previews-container"]}>
          <span className={styles.status}>
            <span className={styles.colour}>{publications.length}</span>
            &nbsp;
            {publications.length === 1 ? "publication" : "publications"}
            &nbsp;found for the below variant
          </span>
          <div className={styles["variant-summary"]}>
            <VariantSummary variant={props.variant} />
          </div>
          <div className={styles["publication-previews"]}>
            {!props.isPhenotypePublicationPage && publications.length > 0 && (
              <div
                className={styles["publication-archive-link"]}
                onClick={getPublicationUpdates}
              >
                Click here to view Publication updates
              </div>
            )}
            <div className={styles.header}>
              <span className={styles["header-title"]}>Publication Titles</span>
              <DebouncedInput
                onChange={(val) => setFilterValue(val.toString())}
                placeholder={`Search publications...`}
                type="text"
                value={filterValue}
                className={styles.input}
              />
            </div>

            <div className={styles["publication-preview-contents"]}>
              {(filterValue.length > 0
                ? publications.filter(
                    (p) =>
                      containsFilterValue(p.title) ||
                      containsFilterValue(p.doi) ||
                      containsFilterValue(p.abstract) ||
                      (p.authors && containsFilterValue(p.authors.join(" "))) ||
                      containsFilterValue(p.journal) ||
                      containsFilterValue(p.pmid?.toString())
                  )
                : publications
              ).map((publication) => (
                <PublicationPreview
                  data={publication}
                  onRemoveClickCallback={() =>
                    setPublicationToBeRemoved(publication)
                  }
                />
              ))}
            </div>
          </div>
        </div>
      )}

      {showPublicationArchiveModal && (
        <Modal
          title="Publication updates archive"
          isClosable={true}
          modalContainerStyle={styles.modal}
          onCloseIconClickCallback={() => setShowPublicationArchiveModal(false)}
        >
          {publicationUpdates.length > 0 ? (
            <div className={styles["publication-updates-container"]}>
              <div className={styles["publication-updates-description"]}>
                <p>
                  Publication updates were checked on every highlighted date
                  shown below.
                </p>
                <p>
                  Dates in bold indicate that new publications were found for
                  this variant. Further details about when publications were
                  added can be found below the calendar.
                </p>
              </div>
              <div className={styles.calendar}>
                <CalendarDisplay
                  markedDates={Array.from(
                    publicationUpdates.map((c) => {
                      const date = c.lastEval.split(" ")[0];

                      return {
                        date: date,
                        update:
                          datesWithUpdates?.find((d) => d === date) ?? false,
                      };
                    })
                  )}
                />
              </div>
              {publicationUpdates
                .filter(
                  (u) =>
                    u.publicationUpdates !== null &&
                    u.publicationUpdates !== undefined
                )
                .map((u) => (
                  <div
                    className={styles["publication-update-info-with-update"]}
                  >
                    <p className={styles["publication-date-checked-container"]}>
                      <span className={styles.bullet}>{"\u25CF"}</span>
                      <span className={styles["publication-date-checked"]}>
                        {u.lastEval}
                      </span>
                    </p>
                    <div className={styles["publication-updates"]}>
                      {u.publicationUpdates.map((p) => (
                        <div
                          className={`${styles["publication-update"]} ${
                            p.link ? styles.link : ""
                          }`}
                          onClick={() => p.link && openInNewWindow(p.link)}
                        >
                          <div className={styles.pub}>
                            <p>
                              <b>Title:</b> {p.title}
                            </p>
                            {p.doi && (
                              <p>
                                <b>DOI:</b> {p.doi}
                              </p>
                            )}
                          </div>
                          <div className={styles["manual-add"]}>
                            {p.isManuallyAdded && (
                              <Icon
                                name="profile"
                                className={styles["pub-archive-profile-icon"]}
                              />
                            )}
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                ))}
            </div>
          ) : (
            <Loader />
          )}
        </Modal>
      )}

      {showAddPublicationModal && (
        <Modal
          title="Add publications"
          isClosable={!isAddingPublications}
          modalContainerStyle={styles.modal}
          onCloseIconClickCallback={() => {
            setShowAddPublicationModal(false);
            setAddedUrls([]);
          }}
        >
          <div className={styles["add-pub-modal-content"]}>
            <p>
              Insert the URLs of the publications you would like to add for this
              variant
            </p>
            <AddPublications
              isAddingPublications={isAddingPublications}
              onPublicationsUpdateCallback={(urls) => setAddedUrls(urls)}
            />
            <Button
              className={styles["add-publications-btn"]}
              disabled={isAddingPublications || addedUrls.length === 0}
              text="Add publication/s"
              onClick={() => addPublications()}
            />
          </div>
          {isAddingPublications && <Loader />}
        </Modal>
      )}

      {publicationToBeRemoved && (
        <Modal modalContainerStyle={styles["confirm-delete-modal"]}>
          <div className={styles["confirm-delete-modal-content"]}>
            <p>
              Are you sure you want to remove Publication&nbsp;
              <b>
                "{publicationToBeRemoved.title}" (DOI:{" "}
                {publicationToBeRemoved.doi})
              </b>
              &nbsp;from this variant ?
            </p>
            <div className={styles["option-btns"]}>
              <Button
                text="Yes, remove it"
                onClick={removePublication}
                className={styles["delete-btn"]}
                disabled={isRemovingPublication}
              />
              <Button
                text="No, return to publication page"
                onClick={closeRemovePublicationModal}
                disabled={isRemovingPublication}
              />
            </div>
          </div>
          {isRemovingPublication && <Loader />}
        </Modal>
      )}
    </div>
  );

  function containsFilterValue(val: string): boolean {
    if (val && val.length > 0) {
      return val?.toLowerCase().includes(filterValue?.toLowerCase());
    }
    return false;
  }

  function getPublicationUpdates() {
    setShowPublicationArchiveModal(true);

    if (publicationUpdates.length === 0) {
      props.vusService
        .getPublicationUpdates({ variantId: parseInt(props.variantId) })
        .then((res) => {
          if (res.isSuccess) {
            setPublicationUpdates(res.variantPublicationUpdates);
            if (res.datesWithUpdates) setDatesWithUpdates(res.datesWithUpdates);
          }
        });
    }
  }

  function addPublications() {
    if (addedUrls.length > 0) {
      setIsAddingPublications(true);
      props.vusService
        .addPublications({
          variantId: props.variantId,
          publicationUrls: addedUrls,
        })
        .then((res) => {
          if (res.isSuccess) {
            res.publications &&
              setPublications(convertPubDates(res.publications));
            setShowAddPublicationModal(false);
            setAddedUrls([]);
          }
          setIsAddingPublications(false);
        });
    }
  }

  function closeRemovePublicationModal() {
    setPublicationToBeRemoved(undefined);
    setIsRemovingPublication(false);
  }

  function removePublication() {
    setIsRemovingPublication(true);

    vusService
      .removePublication({
        variantId: props.variantId,
        publicationId: publicationToBeRemoved.publicationId,
      })
      .then((res) => {
        if (res.isSuccess) {
          res.publications &&
            setPublications(convertPubDates(res.publications));

          closeRemovePublicationModal();
        }
      });
  }
};

export default PublicationViewPage;
