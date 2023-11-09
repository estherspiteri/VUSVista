import React, { useEffect, useState } from "react";
import styles from "./publication-search.module.scss";
import PublicationPreview from "../publication-preview/publication-preview";
import {
  IPublicationPreview,
  IPublicationSearch,
} from "../../models/publication-search/publication-search.model";
import { PublicationService } from "../../services/publication/publication.service";

type PublicationSearchProps = { publicationService?: PublicationService };

const PublicationSearch: React.FunctionComponent<PublicationSearchProps> = (
  props: PublicationSearchProps
) => {
  const [data, setData] = useState<IPublicationSearch | undefined>(undefined);
  const [rsid, setRsid] = useState("");
  const [isSearching, setIsSearching] = useState(false);

  useEffect(() => {
    if (isSearching) {
      props.publicationService?.getPublications({ rsid }).then((res) => {
        if (res.isSuccess) {
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
      <div className={styles["search-container"]}>
        <input
          placeholder="Type variant RSID"
          value={rsid}
          onChange={(e) => setRsid(e.currentTarget.value)}
          disabled={isSearching}
        />
        <div
          className={styles.icon}
          onClick={() =>
            !isSearching && rsid.length > 0 && setIsSearching(true)
          }
        >
          <svg viewBox="0 0 24 24" fill="none">
            <path
              d="M11 6C13.7614 6 16 8.23858 16 11M16.6588 16.6549L21 21M19 11C19 15.4183 15.4183 19 11 19C6.58172 19 3 15.4183 3 11C3 6.58172 6.58172 3 11 3C15.4183 3 19 6.58172 19 11Z"
              stroke-width="2"
              stroke-linecap="round"
              stroke-linejoin="round"
            />
          </svg>
        </div>
      </div>
      {isSearching ? (
        <div className={styles["lds-dual-ring"]}></div>
      ) : (
        <>
          {data && (
            <>
              {data?.publications && (
                <div className={styles["publication-previews"]}>
                  {/* <div className={styles["table-header"]}>
              <span className={styles.pmid}>PMID</span>
              <span>Title</span>
            </div> */}

                  <div className={styles["publication-preview-contents"]}>
                    {data.publications.map((publication) => (
                      <PublicationPreview data={publication} />
                    ))}
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
