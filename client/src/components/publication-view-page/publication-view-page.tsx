import React, { useState } from "react";
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
import { VusService } from "../../services/vus/vus.service";
import { openInNewWindow } from "../../helpers/open-links";
import Button from "../../atoms/button/button";
import Text from "../../atoms/text/text";
import Icon from "../../atoms/icons/icon";
import { convertPubDates } from "../../helpers/date-helper";

type PublicationViewPageProps = {
  description?: string;
  variantId: string;
  variant: IVUSSummary;
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
  const [url, setUrl] = useState("");
  const [addedUrls, setAddedUrls] = useState<string[]>([]);
  const [isAddingPublications, setIsAddingPublications] = useState(false);

  return (
    <div className={styles["publication-view-container"]}>
      <div className={styles["title-container"]}>
        <div className={styles["title-section"]}>
          <div className={styles.title}>Publications</div>
          <Button
            text="Add Publications"
            icon="add"
            onClick={() => setShowAddPublicationModal(true)}
          />
        </div>

        {props.description && (
          <div
            className={styles.description}
            dangerouslySetInnerHTML={{ __html: props.description }}
          />
        )}
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
            <div
              className={styles["publication-archive-link"]}
              onClick={getPublicationUpdates}
            >
              Click here to view Publication updates
            </div>
            <div className={styles.header}>
              <span className={styles["header-title"]}>Publication Titles</span>
              <DebouncedInput
                onChange={(val) => setFilterValue(val.toString())}
                placeholder={`Search publication titles...`}
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
                      containsFilterValue(p.pmid)
                  )
                : publications
              ).map((publication) => (
                <PublicationPreview data={publication} />
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
                          datesWithUpdates?.find((d) => d == date) ?? false,
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
                            {p.isManuallyAdded && <Icon name="profile" />}
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
            <div className={styles["urls-container"]}>
              <div className={styles["new-input"]}>
                <div className={styles.text}>
                  <Text
                    disabled={isAddingPublications}
                    placeholder="Insert URL ..."
                    value={url}
                    autoFocus={true}
                    onChange={(e) => setUrl(e.currentTarget.value)}
                  />
                </div>
                <div className={styles["icon-wrapper"]}>
                  <Icon
                    name="add"
                    className={styles.add}
                    onClick={() => {
                      if (!isAddingPublications && url.trim().length > 0) {
                        if (!addedUrls.find((u) => u === url)) {
                          setAddedUrls(addedUrls.concat(url));
                        }
                        setUrl("");
                      }
                    }}
                  />
                </div>
              </div>
              {addedUrls.map((u) => (
                <div className={styles["added-url"]}>
                  <a href={u}>{u}</a>
                  <Icon
                    name="close"
                    stroke="#008080"
                    onClick={() => {
                      setAddedUrls(addedUrls.filter((url) => url !== u));
                    }}
                  />
                </div>
              ))}
            </div>
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
};

export default PublicationViewPage;
