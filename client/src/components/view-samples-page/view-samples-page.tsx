import React, { useEffect, useState } from "react";
import styles from "./view-samples-page.module.scss";
import Loader from "../../atoms/loader/loader";
import SampleTable from "./sample-table/sample-table";
import { ISample } from "../../models/view-samples.model";
import { SampleService } from "../../services/sample/sample.service";

type ViewAllSamplesProps = { sampleService: SampleService };

const ViewAllSamples: React.FunctionComponent<ViewAllSamplesProps> = (
  props: ViewAllSamplesProps
) => {
  const [sampleList, setSampleList] = useState<ISample[]>(undefined);
  const [selectedSample, setSelectedSample] = useState<ISample>(undefined);

  useEffect(() => {
    props.sampleService?.loadAllSamples().then((res) => {
      if (res.isSuccess) {
        setSampleList(res.sampleList);
      } else {
        //TODO: Handle error
      }
    });
  }, []);
  console.log(typeof selectedSample?.dateOfFileUpload);
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
              setSelectedSample(
                sampleList.filter((s) => s.sampleId === sampleId)[0]
              )
            }
          />
          <div className={styles["sample-info"]}>
            <p className={styles["sample-info-title"]}>Sample information</p>
            {selectedSample && (
              <>
                <div className={styles.info}>
                  <div className={styles["info-title"]}>Genome Version:</div>
                  {selectedSample.genomeVersion}
                </div>

                <div className={styles.info}>
                  <div className={styles["info-title"]}>File upload name:</div>
                  {selectedSample.fileUploadName}
                </div>
                <div className={styles.info}>
                  <div className={styles["info-title"]}>Date uploaded:</div>
                  <span>{`${selectedSample.dateOfFileUpload.getDate()}/${
                    selectedSample.dateOfFileUpload.getMonth() + 1
                  }/${selectedSample.dateOfFileUpload.getFullYear()}`}</span>
                </div>
                {/* <div
                  className={`${styles.info} ${styles["variant-description"]}`}
                >
                  <div className={styles["info-title"]}>Description:</div>
                  {selectedSample.description}
                </div> */}
                <div className={`${styles.info} ${styles.variants}`}>
                  <div className={styles["info-title"]}>Variants:</div>
                  <div className={styles["variant-info"]}>
                    <div className={styles["variant-titles"]}>
                      <div className={styles["variant-id"]}>Variant Ids</div>
                      <div>Genotypes</div>
                      {/* <div>Assigned ACMG rules</div> */}
                    </div>
                    {selectedSample.variants.map((v) => {
                      return (
                        <div className={styles.variant}>
                          <div className={styles["variant-id"]}>
                            {v.variantId}
                          </div>
                          <div>{v.genotype}</div>
                          {/* <div className={styles.rules}>
                        <div className={styles["pathogenic-supporting"]}>
                          PP1
                        </div>
                        <div className={styles["pathogenic-medium"]}>PM6</div>
                        <div className={styles["pathogenic-strong"]}>PS2</div>
                      </div> */}
                        </div>
                      );
                    })}
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
