import React, { useEffect, useState } from "react";
import styles from "./view-samples-page.module.scss";
import { IVus } from "../../models/view-vus.model.tsx/view-vus.model";
import { VusService } from "../../services/vus/vus.service";
import Loader from "../loader/loader";
import SampleTable from "./sample-table/sample-table";

type ViewAllSamplesProps = { vusService: VusService };

const ViewAllSamples: React.FunctionComponent<ViewAllSamplesProps> = (
  props: ViewAllSamplesProps
) => {
  const [vusList, setVusList] = useState<IVus[]>(undefined);
  const [selectedSample, setSelectedSample] = useState<{
    id: string;
    dateCollected: string;
    description: string;
    genomeVersion: string;
    // variants: IVus[];
    fileName: string;
    numOfVariants: number;
  }>(undefined);

  const sampleList = [
    {
      id: "1",
      dateCollected: "22/02/2022",
      description:
        "Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industry's standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book.",
      genomeVersion: "GRCh37",
      fileName: "fileName.xlsx",
      numOfVariants: 3,
    },
  ];

  useEffect(() => {
    props.vusService?.loadAllVus().then((res) => {
      if (res.isSuccess) {
        setVusList(res.vusList);
      } else {
        //TODO: Handle error
      }
    });
  }, []);

  return (
    <div className={styles["view-all-samples-container"]}>
      <div className={styles.title}>Sample List</div>
      <div className={styles.description}>
        <p>
          Below you can find a list of all the samples stored within our
          database.
        </p>
      </div>
      {sampleList ? (
        <div className={styles["samples-container"]}>
          <SampleTable
            sampleList={sampleList}
            onSampleClickCallback={(sampleId) =>
              setSelectedSample(sampleList.filter((s) => s.id === sampleId)[0])
            }
          />
          <div className={styles["sample-info"]}>
            <p className={styles["sample-info-title"]}>Sample information</p>
            {selectedSample && (
              <>
                <div className={styles.info}>
                  <div className={styles["info-title"]}>Genome Version:</div>
                  GRCh37
                </div>

                <div className={styles.info}>
                  <div className={styles["info-title"]}>File upload name:</div>
                  test-samples.xlsx
                </div>
                <div className={styles.info}>
                  <div className={styles["info-title"]}>Date uploaded:</div>
                  02/03/2022
                </div>
                <div
                  className={`${styles.info} ${styles["variant-description"]}`}
                >
                  <div className={styles["info-title"]}>Description:</div>
                  {selectedSample.description}
                </div>
                <div className={`${styles.info} ${styles.variants}`}>
                  <div className={styles["info-title"]}>Variants:</div>
                  <div className={styles["variant-info"]}>
                    <div className={styles["variant-titles"]}>
                      <div className={styles["variant-id"]}>Variant Ids</div>
                      <div>Assigned ACMG rules</div>
                    </div>
                    <div className={styles.variant}>
                      <div className={styles["variant-id"]}>1</div>
                      <div className={styles.rules}>
                        <div className={styles["pathogenic-supporting"]}>
                          PP1
                        </div>
                        <div className={styles["pathogenic-medium"]}>PM6</div>
                        <div className={styles["pathogenic-strong"]}>PS2</div>
                      </div>
                    </div>
                    <div className={styles.variant}>
                      <div className={styles["variant-id"]}>32</div>
                      <div className={styles.rules}>
                        <div className={styles["benign-strong"]}>BS4</div>
                        <div className={styles["benign-supporting"]}>BP6</div>
                      </div>
                    </div>
                    <div className={styles.variant}>
                      <div className={styles["variant-id"]}>24</div>
                      <div className={styles.rules}></div>
                    </div>
                    <div className={styles.variant}>
                      <div className={styles["variant-id"]}>15</div>
                      <div className={styles.rules}>
                        <div className={styles["benign-supporting"]}>BP6</div>
                        <div className={styles["pathogenic-medium"]}>PM6</div>
                      </div>
                    </div>
                  </div>
                </div>
              </>
            )}
          </div>
        </div>
      ) : (
        <Loader />
      )}
    </div>
  );
};

export default ViewAllSamples;
