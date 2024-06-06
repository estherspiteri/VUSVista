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
  const [showPublicationArchiveModal, setShowPublicationArchiveModal] =
    useState(false);
  const [publicationUpdates, setPublicationUpdates] = useState<
    IVariantPublicationUpdates[]
  >([]);

  return (
    <div className={styles["publication-view-container"]}>
      <div className={styles["title-container"]}>
        <div className={styles.title}>Publications</div>
        {props.description && (
          <div
            className={styles.description}
            dangerouslySetInnerHTML={{ __html: props.description }}
          />
        )}
      </div>

      {props.publications && (
        <>
          {props.publications && (
            <div className={styles["publications-previews-container"]}>
              <span className={styles.status}>
                <span className={styles.colour}>
                  {props.publications.length}
                </span>
                &nbsp;
                {props.publications.length === 1
                  ? "publication"
                  : "publications"}
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
                  <span className={styles["header-title"]}>
                    Publication Titles
                  </span>
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
                    ? props.publications.filter(
                        (p) =>
                          containsFilterValue(p.title) ||
                          containsFilterValue(p.doi) ||
                          containsFilterValue(p.abstract) ||
                          (p.authors &&
                            containsFilterValue(p.authors.join(" "))) ||
                          containsFilterValue(p.journal) ||
                          containsFilterValue(p.pmid)
                      )
                    : props.publications
                  ).map((publication) => (
                    <PublicationPreview data={publication} />
                  ))}
                </div>
              </div>
            </div>
          )}
        </>
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
              <p>
                Publication updates were checked on every highlighted date shown
                below.
              </p>
              <p>
                Dates in bold indicate that new publications were found for this
                variant. Further details about when publications were added can
                be found below the calendar.
              </p>
              <div className={styles.calendar}>
                <CalendarDisplay
                  markedDates={Array.from(
                    publicationUpdates.map((c) => {
                      return {
                        date: c.lastEval.split(" ")[0],
                        update:
                          c.publicationUpdates !== undefined &&
                          c.publicationUpdates !== null,
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
                          <p>
                            <b>Title:</b> {p.title}
                          </p>
                          {p.doi && (
                            <p>
                              <b>DOI:</b> {p.doi}
                            </p>
                          )}
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
          }
        });
    }
  }
};

export default PublicationViewPage;
