import React, { useEffect, useState } from "react";
import styles from "./publication-search-page.module.scss";
import PublicationPreview from "./publication-preview/publication-preview";
import {
  IPublicationPreview,
  IPublicationSearch,
} from "../../models/publication-search.model";
import { PublicationService } from "../../services/publication/publication.service";
import Loader from "../../atoms/loader/loader";
import Icon from "../../atoms/icon/icon";

type PublicationSearchProps = { publicationService?: PublicationService };

const PublicationSearch: React.FunctionComponent<PublicationSearchProps> = (
  props: PublicationSearchProps
) => {
  const [data, setData] = useState<IPublicationSearch | undefined>(undefined);
  const [rsid, setRsid] = useState("");
  const [isSearching, setIsSearching] = useState(false);
  const [searchedRsid, setSearchedRsid] = useState("");

  //read rsid from query parameter and start search
  useEffect(() => {
    const urlParams = new URLSearchParams(window.location.search);
    const rsid = urlParams.get("rsid");

    if (rsid?.length > 0) {
      setRsid(rsid);
      setIsSearching(true);
    }
  }, []);

  useEffect(() => {
    if (isSearching) {
      props.publicationService?.getPublications({ rsid }).then((res) => {
        if (res.isSuccess) {
          setSearchedRsid(rsid);
          setData({
            ...res.publicationSearch,
            publications: convertPubDates(res.publicationSearch.publications),
          });
          setIsSearching(false);
        } else {
          //TODO: handle error
        }
      });
    }
  }, [isSearching, rsid]);

  return (
    <div className={styles["publication-search-container"]}>
      <div className={styles["title-container"]}>
        <div className={styles.title}>Publications Search</div>
        <div className={styles.description}>
          <p>Look up publications for a particular VUS using its RSID.</p>
        </div>
      </div>
      <div className={styles["search-container"]}>
        <input
          autoFocus={true}
          placeholder="Type variant RSID"
          value={rsid}
          onChange={(e) => setRsid(e.currentTarget.value)}
          disabled={isSearching}
        />
        {/** TODO: add onEnter = click */}
        <div
          className={`${styles.icon} ${
            isSearching ? styles["searching-icon"] : ""
          }`}
          onClick={() =>
            !isSearching && rsid.length > 0 && setIsSearching(true)
          }
        >
          <Icon name="search" />
        </div>
      </div>
      {isSearching ? (
        <Loader />
      ) : (
        <>
          {data && (
            <>
              {data?.publications && (
                <div className={styles["publications-previews-container"]}>
                  <span className={styles.status}>
                    <span className={styles.colour}>
                      {data.publications.length}
                    </span>
                    &nbsp;
                    {data.publications.length === 1
                      ? "publication"
                      : "publications"}
                    &nbsp;found for&nbsp;
                    <span className={styles.colour}>{searchedRsid}</span>
                  </span>
                  <div className={styles["publication-previews"]}>
                    <div className={styles.header}>
                      <span className={styles.pmid}>PMID</span>
                      <span>Title</span>
                    </div>

                    <div className={styles["publication-preview-contents"]}>
                      {data.publications.map((publication) => (
                        <PublicationPreview data={publication} />
                      ))}
                    </div>
                  </div>
                </div>
              )}

              {!data.isLitvarIdFound && <h1>Litvar Id not found!!</h1>}
            </>
          )}
        </>
      )}
    </div>
  );

  function convertPubDates(
    publications: IPublicationPreview[]
  ): IPublicationPreview[] {
    let updatedPublications: IPublicationPreview[] = [];

    if (publications) {
      //format publication date
      publications.forEach((publication) => {
        let pub = publication;
        pub.date = new Date(publication.date);
        updatedPublications = updatedPublications.concat(pub);
      });
    }

    return updatedPublications;
  }
};

export default PublicationSearch;
