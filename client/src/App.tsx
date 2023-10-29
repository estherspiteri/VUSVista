import React, { useEffect, useState } from "react";
import styles from "./App.module.scss";
import SearchSVG from "./search-circle-svgrepo-com.svg";
import LiteraturePreview from "./components/literature-preview/literature-preview";

type AppProps = {};

const App: React.FunctionComponent<AppProps> = () => {
  const [data, setData] = useState(undefined);
  const [rsid, setRsid] = useState("");
  const [isSearching, setIsSearching] = useState(false);

  useEffect(() => {
    if (isSearching) {
      fetch(`/litvar/publications/${rsid}`, {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      })
        .then((response) => {
          if (response.status === 200) return response.text();
          //TODO: else show error page
        })
        .then((res) => {
          return JSON.parse(res);
        })
        .then((data) => {
          setData(Object.keys(data.data).length == 0 ? undefined : data.data);
          setIsSearching(false);
        })
        .catch((error) => console.error("error============:", error));
    }
  }, [isSearching, rsid]);

  return (
    <div className={styles["container"]}>
      <div className={styles["literature-search-container"]}>
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
              <div className={styles["literature-previews"]}>
                {/* <div className={styles["table-header"]}>
                  <span className={styles.pmid}>PMID</span>
                  <span>Title</span>
                </div> */}

                <div className={styles["literature-preview-contents"]}>
                  {Object.values(data?.publications).map((publication) => (
                    <LiteraturePreview
                      pmid={publication.pmid}
                      title={publication.title}
                      date={new Date(publication.date)}
                      authors={publication.authors}
                      journal={publication.journal}
                      abstract={publication.abstract}
                      isSupplementaryMaterialMatch={
                        publication.is_sup_mat_match
                      }
                    />
                  ))}
                </div>
              </div>
            )}
          </>
        )}
      </div>
    </div>
  );
};

export default App;
